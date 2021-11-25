#!/usr/bin/env python
## Simply merge the overlapping regions. Using a more precise boundary to decide the break point
## This version is used for results from 'MetaSV_v5' version.

import os
import sys
import math
import re

'''
defind functions
========================================================================================================================
'''
def sortByPos(x,y,key=1):
	'''
	sort variants_chr[chr] by position (key=1), asending
	'''
	if int(x[key])<int(y[key]):
		return -1
	elif int(x[key])>int(y[key]):
		return 1
	else:
		return 0

def modifyLine(var,pos,end,ref,filter,prec,rd_only):
	line = var.split('\t')
	line[INFO_col] = line[INFO_col].replace('IMPRECISE=.;','')
	imprecise = '' if prec else ';IMPRECISE=.'
	info_match = INFOp.search(line[INFO_col])
	if not info_match:
		sys.exit('1.Check ' + var)
	line_info = info_match.group('i0') + str(end) + imprecise + info_match.group('i1') + info_match.group('svtool') + info_match.group('i2') + info_match.group('sources') + info_match.group('i3') + info_match.group('num') + info_match.group('i4') + info_match.group('meth') + info_match.group('i5')
	new_line = '\t'.join([line[0],str(pos),line[2],ref,line[4],line[5],filter,line_info,line[8],line[9]])
	return [new_line,pos,end,ref,filter,prec,rd_only]

def updateLine(var1, var2, pos, end, ref, filter, prec, rd_only):
	line1 = var1.split('\t')
	line2 = var2.split('\t')
	line1[INFO_col] = line1[INFO_col].replace('IMPRECISE=.;','')
	line2[INFO_col] = line2[INFO_col].replace('IMPRECISE=.;','')
	info1_match = INFOp.search(line1[INFO_col])
	info2_match = INFOp.search(line2[INFO_col])
	if (not info1_match) or (not info2_match):
		sys.exit('2.Check ' + '\t'.join(line1) + ' or ' + '\t'.join(line2))
	imprecise = '' if prec else ';IMPRECISE=.'
	line_sources = info1_match.group('sources') + ',' + info2_match.group('sources')
	two_meth = info1_match.group('meth').split(',')
	two_meth.extend(info2_match.group('meth').split(','))
	line_meth = []
	[line_meth.append(i) for i in two_meth if not i in line_meth]
	num = len(line_meth)
	line_info = info1_match.group('i0') + str(end) + imprecise + info1_match.group('i1') + info1_match.group('svtool') + info1_match.group('i2') + line_sources + info1_match.group('i3') + str(num) + info1_match.group('i4') + ','.join(line_meth) + info1_match.group('i5') + info2_match.group('i5')
	line = '\t'.join([line1[0],str(pos),line1[2],ref,line1[4],line1[5],filter,line_info,line1[8],line1[9]])
	return [line,pos,end,ref,filter,prec,rd_only]

'''
========================================================================================================================
'''

INFOp=re.compile(r'(?P<i0>^.*)(?<=\bEND=)(?P<end>\d+)(?=;)(?P<i1>.*SVTOOL=)(?P<svtool>\w+)(?P<i2>.*SOURCES=)(?P<sources>[0-9A-Za-z\-,]+)(?P<i3>.*NUM_SVMETHODS=)(?P<num>\d+)(?P<i4>.*SVMETHOD=)(?P<meth>[\w,]+)(?P<i5>.*$)')

if len(sys.argv) != 3:
	sys.exit("Argument: MetaSV.vcf Type(DEL|DUP)\nVery Important: Only one SV type!!!")

MetaSV = sys.argv[1]
svtype = sys.argv[2]
WD = os.path.dirname(os.path.abspath(MetaSV))

RD_overlap_cri = 300                    # if only MetaSV-RD supported, the criteria for the overlapping region should > 300bp

## MetaSV Info
POS_col = 2-1
REF_col = 4-1
FILTER_col = 7-1
INFO_col = 8-1

METASV = file(MetaSV,'r')
annotations = []
metasv_header = ''
variants_chr = {}
chrom = []

while True:
	line=METASV.readline().rstrip()
	if len(line) <= 0:
		break
	if line.startswith('##'):
		annotations.append(line + '\n')
	elif line.startswith('#'):
		try:
			if len(metasv_header) > 0:
				raise ExistedHeader
			else:
				metasv_header = line + '\n'
		except ExistedHeader:
			sys.exit('Multiple header in ' + MetaSV + '.')
	else:
		chr=line.partition('\t')[0]
		if chr.upper() == 'X':
			chr='23'
			line = '23' + '\t' + line.partition('\t')[2]
		elif chr.upper() == 'Y':

			chr='24'
			line = '24' + '\t' + line.partition('\t')[2]
		filter = line.split('\t')[FILTER_col]
		info = line.split('\t')[INFO_col]
		pos = line.split('\t')[POS_col]
		ref = line.split('\t')[REF_col]
		end = str()
#		if info.startswith('END='):
#			end = info.partition(';')[0].partition('=')[2]
#		if INFOp.search(line).group('end'):		# modified by Fu Ruiqing, 2016.11.18
		if INFOp.search(line):
			end = INFOp.search(line).group('end')
		else:
			if info.find(';SVTYPE=INV;') != -1 or info.find(';SVTYPE=INS;') != -1:
				continue
			else:
				# like ITX
				sys.exit("No END info ,check " + line)
		PRECISE = True if info.find('IMPRECISE')==-1 else False
		RD_only = True if info.find(';SVMETHOD=RD;')!=-1 else False
		if not chr in chrom:
			chrom.append(chr)
			variants_chr[chr] = [[line,pos,end,ref,filter,PRECISE,RD_only]]
		else:
			if int(pos) < int(variants_chr[chr][-1][1]):
				sys.exit(MetaSV + ' should be sorted first!')
			else:
				variants_chr[chr].append([line,pos,end,ref,filter,PRECISE,RD_only])
METASV.close()

output = []
for chr in chrom:
	chr_sv = variants_chr[chr][:]
	while chr_sv:
		curSV = chr_sv.pop(0)
		tmp_pos = int(curSV[1])
		tmp_end = int(curSV[2])
		tmp_ref = curSV[3]
		tmp_filter = curSV[4]
		tmp_prec = curSV[5]
		tmp_rd = curSV[6]
		try:
			while int(chr_sv[0][1]) <= tmp_end:
				if (tmp_end - int(chr_sv[0][1]) + 1 <= RD_overlap_cri) and tmp_rd:
					if int(chr_sv[0][1])-1 - tmp_pos + 1>=100:
						tmp_end = int(chr_sv[0][1])-1
						curSV = modifyLine(curSV[0],tmp_pos,tmp_end,tmp_ref,tmp_filter,tmp_prec,tmp_rd)
					else:
						curSV = []
					break	# break while
				elif (tmp_end - int(chr_sv[0][1]) + 1 <= RD_overlap_cri) and chr_sv[0][6]:
					if int(chr_sv[0][2])-(tmp_end+1)+1 >= 100:
						chr_sv[0]=modifyLine(chr_sv[0][0],tmp_end+1,chr_sv[0][2],'.',chr_sv[0][4],chr_sv[0][5],chr_sv[0][6])
						chr_sv.sort(sortByPos)
					else:
						rm = chr_sv.pop(0)
						del rm
				else:
					nextSV = chr_sv.pop(0)
					## priority(PRECISE) > priority(PASS)
					if tmp_prec and nextSV[5]:
						tmp_pos = min(tmp_pos,int(nextSV[1]))
						tmp_end = max(tmp_end,int(nextSV[2]))
						tmp_filter = 'PASS'
						tmp_ref = tmp_ref if tmp_pos <= int(nextSV[1]) else nextSV[3]
						tmp_rd = tmp_rd if not tmp_rd else nextSV[6]
					elif tmp_prec:
						tmp_pos = tmp_pos
						tmp_end = tmp_end
						tmp_ref = tmp_ref
						tmp_filter = 'PASS'
						tmp_prec = tmp_prec
						tmp_rd = tmp_rd if not tmp_rd else nextSV[6]
					elif nextSV[5]:
						tmp_pos = int(nextSV[1])
						tmp_end = int(nextSV[2])
						tmp_ref = nextSV[3]
						tmp_filter = 'PASS'
						tmp_prec = nextSV[5]
						tmp_rd = tmp_rd if not tmp_rd else nextSV[6]
					elif tmp_filter == 'PASS' and nextSV[4] == 'PASS':
						tmp_pos = min(tmp_pos,int(nextSV[1]))
						tmp_end = max(tmp_end,int(nextSV[2]))
						tmp_filter = 'PASS'
						tmp_ref = tmp_ref if tmp_pos <= int(nextSV[1]) else nextSV[3]
						tmp_rd = tmp_rd if not tmp_rd else nextSV[6]
					elif tmp_filter == 'PASS':
						tmp_pos = tmp_pos
						tmp_end = tmp_end
						tmp_ref = tmp_ref
						tmp_filter = 'PASS'
						tmp_prec = tmp_prec
						tmp_rd = tmp_rd if not tmp_rd else nextSV[6]
					elif nextSV[4] == 'PASS':
						tmp_pos = int(nextSV[1])
						tmp_end = int(nextSV[2])
						tmp_ref = nextSV[3]
						tmp_filter = 'PASS'
						tmp_prec = nextSV[5]
						tmp_rd = tmp_rd if not tmp_rd else nextSV[6]
					else:
						tmp_pos = min(tmp_pos,int(nextSV[1]))
						tmp_end = max(tmp_end,int(nextSV[2]))
						tmp_ref = tmp_ref if tmp_pos <= int(nextSV[1]) else nextSV[3]
						tmp_rd = tmp_rd if not tmp_rd else nextSV[6]
					curSV = updateLine(curSV[0],nextSV[0],tmp_pos,tmp_end,tmp_ref,tmp_filter,tmp_prec,tmp_rd)
			if curSV:
				output.append(curSV[0]+'\n')
		except IndexError:
			output.append(curSV[0]+'\n')

noOverlap = WD + '/' + os.path.basename(MetaSV).rpartition('.')[0] + '.nonOverlap.vcf'
OUTPUT=file(noOverlap,'w')
OUTPUT.writelines(annotations)
OUTPUT.writelines(metasv_header)
OUTPUT.writelines(output)
OUTPUT.close()
