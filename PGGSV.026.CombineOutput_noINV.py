#coding: utf-8
#2020-7-23
#该脚本是用于合并metasv_manta结果里的DEL、DUP、INS结果

import sys
import gzip 

ref_DEL_DUP_path = sys.argv[1]
ref_INS_path = sys.argv[2]
fout_SV_path = sys.argv[3]

fout_SV_file = open(fout_SV_path,"w")
ref_tmp_chr_dt = {}
for i in range(1,23):
	chr = str(i)
	ref_tmp_chr_dt[chr] = 0

ref_chr_dt = {}
ref_chr_list = []
with open(ref_DEL_DUP_path) as file:
	for line in file:
		val = line.strip()
		if val.count("#") != 0:
			pass
		elif val.count("#") == 0:
			j = val.split()
			chr = j[0]
			pos = j[1]
			if chr == "23":
				chr = "chrX"
			elif chr == "24":
				chr = "chrY"
			else:
				if chr not in ref_tmp_chr_dt:
					print(ref_DEL_DUP_path)
					break
				chr = "chr" + chr
			j[0] = chr 
			new_line = "\t".join(j) + "\n"
			if chr not in ref_chr_dt:
				ref_chr_list.append(chr)
				ref_chr_dt[chr] = []
				ref_chr_dt[chr].append([int(pos),new_line])
			elif chr in ref_chr_dt:
				ref_chr_dt[chr].append([int(pos),new_line])


with gzip.open(ref_INS_path) as file:
	for line in file:
		val = line.strip()
		if val.count("#") != 0:
			fout_SV_file.write(line)
		else:
			j = val.split()
			chr = j[0]
			pos = j[1]
			SV_type_1 = j[2]
			SV_type_2 = j[4]
			if (SV_type_1.count("INS") != 0) and (SV_type_2.count("INS") == 0):
				if chr in ref_chr_dt:
					ref_chr_dt[chr].append([int(pos),line])

for chr in ref_chr_list:
	fout_sort_chr_list = sorted(ref_chr_dt[chr],key=lambda x:int(x[0]))
	for nu in range(len(fout_sort_chr_list)):
		fout_line = fout_sort_chr_list[nu][1]
		fout_SV_file.write(fout_line)

fout_SV_file.close()
