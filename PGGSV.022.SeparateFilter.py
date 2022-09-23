#!/usr/bin/env python
#####CNVs that are to be kept follows:
##### 1. size_low_cri <= Size <= size_up_cri
##### 2. PASS or (resulted from CNVnator and CNVnator_size_cri(1kb) <= Size <= size_up_cri.

import os
import sys
import re
import gzip

if len(sys.argv) < 5 or len(sys.argv) > 7:
	sys.exit("Argument: CNV.vcf[.gz] BED_file(H) CNVfraction(>) RD_confSize '-'Size_low(>=1bp) '+'Size_up(<=3e9)\nNote1: Get the SVtype from the INFO column\nNote2: Use '-'/'+' to specify lowwer/upper bound for size, no blanks!\nNote3: The BED_file is typically a gap file.")
vcf = sys.argv[1]
bedfile=sys.argv[2]
frac=float(sys.argv[3])
CNVnator_size_cri = float(sys.argv[4])
size_low_cri = 1
size_up_cri = 3e9
#CNVnator_size_cri = 10e3
#CNVnator_size_cri = 1e3			# changed since 2016.08.11
if len(sys.argv)>5:
	for par in sys.argv[5:]:
		if par.startswith('-'):
			size_low_cri = int(par.lstrip('-'))
		elif par.startswith('+'):
			size_up_cri = int(par.lstrip('+'))
		else:
			sys.exit("Argument: CNV.vcf[.gz] BED_file CNVfraction(>) RD_confSize '-'Size_low(>=1bp) '+'Size_up(<=3e9)\nNote1: Get the SVtype from the INFO column\nNote2: Use '-'/'+' to specify lowwer/upper bound for size, no blanks!\n")
WD = os.path.dirname(os.path.abspath(vcf))
CHR_col = 1-1
POS_col = 2-1
FILT_col = 7-1					# if this variant gets a 'PASS', it will be kept, without considerring the size
INFO_col = 8-1
BED_col = (1-1,2-1,3-1)				# the column of bed file, as chr, start, end
ENDp = re.compile(r'(?<=\bEND=)(?P<end>\d+)(?=;)')
SVTYPEp = re.compile(r'(?<=\bSVTYPE=)(?P<svtype>[\w:\-_]+)(?=;)')
SVLENp =re.compile(r'(?<=\bSVLEN=)(?P<svlen>\-?\d+)(?=;|$)')
#CNVnatorp = re.compile(r'\bSOURCES=[\d\-]+\-CNVnator\b')
CNVnatorp = re.compile(r'CNVnator') # xiebo edited

bedregion={}
try:
	BED=file(bedfile,'r')
	bedheader=BED.readline()
	while True:
		line = BED.readline()
		if len(line) <= 0:
			break
		else:
			chrom = line.split("\t")[BED_col[0]].upper().replace('CHR','')
			chrom = chrom.replace("X",'23')
			chrom = chrom.replace("Y",'24')
			start = int(line.split("\t")[BED_col[1]])
			end = int(line.split("\t")[BED_col[2]])
			if bedregion.has_key(chrom):
				if start < bedregion[chrom][-1][1]:
					print line
					sys.exit("You should first sort the regions of " + bedfile)
				bedregion[chrom].append([chrom,start,end])
			else:
				bedregion[chrom]=[[chrom,start,end]]
	BED.close()
except IOError:
	print "No gap file provided! Just ignore!"

VCF=''
if vcf.endswith('.gz'):
	VCF=gzip.open(vcf,'rb')
else:
	VCF=file(vcf,'r')

annotations = []
SVtypes = []
sep_sv = {}
header = ''

while True:
	line = VCF.readline()
	if len(line) <=0:
		break
	else:
		if line.startswith('##'):
			annotations.append(line)
		elif line.startswith('#'):
			try:
				if len(header) > 0:
					raise ExistedHeader
				else:
					header=line
			except ExistedHeader:
				sys.exit('Multiple header in ' + del_vcf + '.')
		else:
			chrom = line.split('\t')[CHR_col].upper().replace('CHR','')
			chrom = chrom.replace('X','23')
			chrom = chrom.replace('Y','24')
			pos = int(line.split('\t')[POS_col])
			end = int()
			svlen = int()
			svtype = ''
			CNVnatorFlag = True if CNVnatorp.search(line.split('\t')[INFO_col]) else False
			try:
				svtype = SVTYPEp.search(line.split('\t')[INFO_col]).groupdict()['svtype']
			except:
				sys.exit("Cannot find SVTYPE info in\n" + line )
			try:
				end = int(ENDp.search(line.split('\t')[INFO_col]).groupdict()['end'])
			except:
				if svtype == 'INV' or svtype =='INS':
				## some INV or INS type don't have END info, so use 'SVLEN' instead
					try:
						svlen = int(SVLENp.search(line.split('\t')[INFO_col]).groupdict()['svlen'])
					except:
						sys.exit("Cannot find END or SVLEN info in\n" + line)
#					try:
#						sep_sv[svtype].append(line)
#					except KeyError:
#						SVtypes.append(svtype)
#						sep_sv[svtype]=[line]
#					except:
#						sys.exit("I don't know")
				else:
					sys.exit("Cannot find END info in\n" + line )
			else:
##				if line.split('\t')[FILT_col] == 'PASS' or ((int(end) - pos + 1 >= size_low_cri) and (int(end)-pos+1 <=size_up_cri)):
				size = end - pos if end > 0 else svlen
				if abs(size) >= size_low_cri and abs(size) <= size_up_cri:
					if svtype == 'DEL' or svtype == 'DUP':
						if not end > 0:
							sys.exit("Cannot find END info in\n" + line)
						bedregion_filt=False
						if bedregion.has_key(chrom):
							for region in bedregion[chrom]:
								if region[1] > end:
									continue
								elif region[2] < pos:
									break
								else:
									overlap_frac = (min(end,region[2]) - max(pos,region[1]) + 1)/float(size)
									if overlap_frac > frac:
										bedregion_filt=True
						if bedregion_filt:
							continue
					if line.split('\t')[FILT_col] == 'PASS':
						try:
							sep_sv[svtype].append(line)
						except KeyError:
							SVtypes.append(svtype)
							sep_sv[svtype]=[line]
						except:
							sys.exit("I don't know")
					else:
						if CNVnatorFlag and abs(size) >= CNVnator_size_cri:
							try:
								sep_sv[svtype].append(line)
							except KeyError:
								SVtypes.append(svtype)
								sep_sv[svtype]=[line]
							except:
								sys.exit("I don't know")
VCF.close()

for svtype in SVtypes:
	sep_vcf = WD + '/' + os.path.basename(vcf).strip('.gz').rpartition('.')[0] + '.' + svtype + '.vcf'
	SEP=file(sep_vcf,'w')
	SEP.writelines(annotations)
	SEP.writelines(header)
	SEP.writelines(sep_sv[svtype])
	SEP.close()

