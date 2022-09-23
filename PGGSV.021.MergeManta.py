#coding: utf-8

import os 
import sys
import gzip

BIN = str(os.path.split(os.path.realpath(__file__))[0]) + "/"

ref_manta_DEL_DUP_path = sys.argv[1]
ref_meta_sv_path = sys.argv[2]
fout_meta_sv_path = sys.argv[3]

#create chr list and dict
ref_chr_list = []
ref_Third_chr_dt = {}
ref_manta_chr_dt = {}
ref_manta_DEL_chr_dt = {}
ref_manta_DUP_chr_dt = {}
ref_metaSV_chr_dt = {}
ref_metaSV_manta_cutoff_50_chr_dt = {}

for i in range(1,23):
    chr = "chr" + str(i)
    ref_chr_list.append(chr)
    ref_Third_chr_dt[chr] = []
    ref_manta_chr_dt[chr] = []
    ref_manta_DEL_chr_dt[chr] = []
    ref_manta_DUP_chr_dt[chr] = []
    ref_manta_chr_dt[chr] = []
    ref_metaSV_chr_dt[chr] = []
    ref_metaSV_manta_cutoff_50_chr_dt[chr] = []

ref_sex_chr_list = ["chrX","chrY"]
for i in ref_sex_chr_list:
    ref_chr_list.append(str(i))
    ref_Third_chr_dt[str(i)] = []
    ref_manta_chr_dt[str(i)] = []
    ref_manta_DEL_chr_dt[str(i)] = []
    ref_manta_DUP_chr_dt[str(i)] = []
    ref_manta_chr_dt[str(i)] = []
    ref_metaSV_chr_dt[str(i)] = []
    ref_metaSV_manta_cutoff_50_chr_dt[str(i)] = []


fout_metaSV_manta_cutoff_50_head_list = []


#该函数是根据提供的cufoff，判断SV之间是否是重叠
def run_loose(ref_list,other_SV_type,other_pos,other_length,cutoff):
    test_nu = 0
    for i in range(len(ref_list)):
        Tujia_line = ref_list[i]
        Tujia_pos = Tujia_line[1]
        Tujia_end = Tujia_line[2]
        Tujia_length = Tujia_line[3]
        Tujia_SV_type = Tujia_line[4]
        Tujia_start = int(Tujia_pos)
        Tujia_end = int(Tujia_pos) + int(Tujia_length)
        other_start = int(other_pos)
        other_end = int(other_pos) + int(other_length)
        if (Tujia_SV_type.count("DEL") !=0) and (other_SV_type.count("DEL") != 0):
            if other_start <= Tujia_start <= other_end:
                if Tujia_end <= other_end:
                    overlap_length = Tujia_end - Tujia_start
                elif Tujia_end >= other_end:
                    overlap_length = other_end - Tujia_start

                if (cutoff*int(other_length) <= overlap_length) and (cutoff*int(Tujia_length) <= overlap_length):
                    test_nu = 1
                    break
            if Tujia_start <= other_start <= Tujia_end:
                if other_end <= Tujia_end:
                    overlap_length = other_end - other_start
                elif other_end > Tujia_end:
                    overlap_length = Tujia_end - other_start

                if (cutoff*int(other_length) <= overlap_length) and (cutoff*int(Tujia_length) <= overlap_length):
                    test_nu = 1
                    break
        if (Tujia_SV_type.count("DUP") !=0) and (other_SV_type.count("DUP") != 0):
            if other_start <= Tujia_start <= other_end:
                if Tujia_end <= other_end:
                    overlap_length = Tujia_end - Tujia_start
                elif Tujia_end >= other_end:
                    overlap_length = other_end - Tujia_start
                if (cutoff*int(other_length) <= overlap_length) and (cutoff*int(Tujia_length) <= overlap_length):
                    test_nu = 1
                    break
            if Tujia_start <= other_start <= Tujia_end:
                if other_end <= Tujia_end:
                    overlap_length = other_end - other_start
                elif other_end > Tujia_end:
                    overlap_length = Tujia_end - other_start
                if (cutoff*int(other_length) <= overlap_length) and (cutoff*int(Tujia_length) <= overlap_length):
                    test_nu = 1
                    break
    return(test_nu,Tujia_SV_type,Tujia_pos,Tujia_end,Tujia_length)


#读取manta的数据
def read_manta_file(ref_manta_DEL_path,ref_manta_chr_dt):
    with gzip.open(ref_manta_DEL_path) as file:
        for line in file:
            if line.count("#") != 0:
                pass
            else:
                j = line.strip().split()
                filter = j[6]
                if filter == "PASS":
                    end = 0
                    genotype = j[-1].split(":")[0]
                    if genotype.count("1") != 0:
                        chr = j[0]
                        pos = j[1]
                        ann_list = j[7].split(";")
                        SV_type = j[2]
                        if chr in ref_Third_chr_dt:
                            if (SV_type.count("DEL") != 0) or (SV_type.count("DUP") != 0)  :
                                for ann in ann_list:
                                    if ann.count("END") != 0:
                                        if ann.split("=")[0] == "END":
                                            end = ann.split("=")[1]
                                length = int(end) - int(pos)
                                if (int(length) < 10000):
                                    length = str(length)
                                    ref_manta_chr_dt[chr].append([chr,int(pos),end,length,SV_type])

    return(ref_manta_chr_dt)

ref_manta_DEL_DUP_chr_dt = read_manta_file(ref_manta_DEL_DUP_path,ref_manta_DEL_chr_dt)


#把manta的DEL和DUP结果写入同一个字典
ref_manta_chr_dt = {}
for chr in ref_chr_list:
    ref_DEL_DUP_list = ref_manta_DEL_DUP_chr_dt[chr]
    ref_manta_chr_dt[chr] = sorted(ref_DEL_DUP_list,key=lambda x:x[1])


#读取metasv的数据，与三代测序做比较，并且在cutoff为50%或是80%的情况下，与manta数据合并，再分别与三代测序数据做比较
def metasv_line_convert(j,ref_ann_list,ref_ann_dt,manta_cutoff_50_SV_type,manta_cutoff_50_pos,manta_cutoff_50_end,manta_cutoff_50_length):
    chr = j[0]
    pos = j[1]
    FILTER = j[6]
    if FILTER == "LowQual":
        FILTER = "PASS"
    INFO = j[7]
    NUM_SVMETHODS = ref_ann_dt["NUM_SVMETHODS"][1].split("=")[1]
    ref_ann_dt["NUM_SVMETHODS"][1] = "NUM_SVMETHODS=" + str(int(NUM_SVMETHODS) + 1)
    NUM_SVTOOLS = ref_ann_dt["NUM_SVTOOLS"][1].split("=")[1]
    ref_ann_dt["NUM_SVTOOLS"][1] = "NUM_SVTOOLS=" + str(int(NUM_SVTOOLS) + 1)
    SOURCES = ref_ann_dt["SOURCES"][1].split("=")[1]
    manta_source = chr + "-" + str(manta_cutoff_50_pos) + "-" + chr + "-" + str(manta_cutoff_50_end) + "-" + str(manta_cutoff_50_length) + "-Manta"
    SOURCES = SOURCES + "," + manta_source
    ref_ann_dt["SOURCES"][1] = "SOURCES=" +  SOURCES
    if (SOURCES.count("BreakSeq") == 0) and (SOURCES.count("Lumpy") == 0):
        pos = str(manta_cutoff_50_pos)
        ref_ann_dt["END"][1] = "END=" + str(manta_cutoff_50_end)
        if manta_cutoff_50_SV_type.count("DEL") != 0:
            ref_ann_dt["SVLEN"][1] = "SVLEN=-" + str(manta_cutoff_50_length)
        elif manta_cutoff_50_SV_type.count("DUP") != 0:
            ref_ann_dt["SVLEN"][1] = "SVLEN=" + str(manta_cutoff_50_length)

    ann_list = j[7].split(";")
    for ann_index in ref_ann_list:
        ann_nu = ref_ann_dt[ann_index][0]
        ann = ref_ann_dt[ann_index][1]
        ann_list[ann_nu] = ann
    INFO = ";".join(ann_list)
    fout_list = [j[0],pos,j[2],j[3],j[4],j[5],FILTER,INFO,j[8],j[9]]
    return(fout_list)


with open(ref_meta_sv_path) as file:
    for line in file:
        val = line.strip()
        if val.count("#") != 0:
            fout_metaSV_manta_cutoff_50_head_list.append(line)
        else:
            ref_ann_dt = {}
            ref_ann_list = ["END","SVLEN","SOURCES","NUM_SVMETHODS","NUM_SVTOOLS"]
            j = val.split()
            chr = j[0]
            pos = int(j[1])
            SV_type = j[4]
            INFO = j[7]
            ann_list = INFO.split(";")
            for nu in range(len(ann_list)):
                ann = ann_list[nu]

                if ann.split("=")[0] == "END":
                    end = int(ann.split("=")[1])
                    length = end - pos
                    ref_ann_dt["END"] = [nu,ann]
                elif ann.split("=")[0] == "SVLEN":
                    ref_ann_dt["SVLEN"] = [nu,ann]
                elif ann.split("=")[0] == "SOURCES":
                    ref_ann_dt["SOURCES"] = [nu,ann]
                elif ann.split("=")[0] == "NUM_SVMETHODS":
                    ref_ann_dt["NUM_SVMETHODS"] = [nu,ann]
                elif ann.split("=")[0] == "NUM_SVTOOLS":
                    ref_ann_dt["NUM_SVTOOLS"] = [nu,ann]
            if (SV_type.count("DEL") == 0) and (SV_type.count("DUP") == 0):
                ref_metaSV_manta_cutoff_50_chr_dt[chr].append(j)
            else:
                if length >= 10000:
                    ref_metaSV_manta_cutoff_50_chr_dt[chr].append(j)
                elif length < 10000:
                    ref_manta_list = ref_manta_chr_dt[chr]
                    if len(ref_manta_list) == 0:
                        ref_metaSV_manta_cutoff_50_chr_dt[chr].append(j)
                    else:
                        metasv_manta_cutoff_50_test_nu,manta_cutoff_50_SV_type,manta_cutoff_50_pos,manta_cutoff_50_end,manta_cutoff_50_length = run_loose(ref_manta_list,SV_type,pos,length,0.5)
                        if metasv_manta_cutoff_50_test_nu == 1:
                            metasv_manta_cutoff_50_list = metasv_line_convert(j,ref_ann_list,ref_ann_dt,manta_cutoff_50_SV_type,manta_cutoff_50_pos,manta_cutoff_50_end,manta_cutoff_50_length)
                            ref_metaSV_manta_cutoff_50_chr_dt[chr].append(metasv_manta_cutoff_50_list)
                        elif metasv_manta_cutoff_50_test_nu == 0:
                            ref_metaSV_manta_cutoff_50_chr_dt[chr].append(j)


fout_meta_manta_cutoff_50_sv_file = open(fout_meta_sv_path,"w")

# fout fout_metaSV_manta_cutoff_50
for fout_line in fout_metaSV_manta_cutoff_50_head_list:
    fout_meta_manta_cutoff_50_sv_file.write(fout_line)


for chr in ref_chr_list:
    fout_metaSV_manta_cutoff_50_list = sorted(ref_metaSV_manta_cutoff_50_chr_dt[chr],key=lambda x:int(x[1]))
    for fout_metaSV_manta_cutoff_50 in fout_metaSV_manta_cutoff_50_list:
        result = tuple(fout_metaSV_manta_cutoff_50)
        STR_result = "\t".join(result)
        fout_meta_manta_cutoff_50_sv_file.write(STR_result)
        fout_meta_manta_cutoff_50_sv_file.write("\n")