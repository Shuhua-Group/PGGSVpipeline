#!/usr/bin/env python

import os
import sys
import re
import gzip

if len(sys.argv) < 3 or len(sys.argv) > 4:
	sys.exit("Argument: CNV.vcf[.gz] Size_low(bp) [Size_up(<=3e9)]\nNote1: Get the SVtype from the INFO column")
vcf = sys.argv[1]
size_low_cri = int(sys.argv[2])
size_up_cri = 3e9
size_suff = 'size_' + str(size_low_cri)
if len(sys.argv)>3:
	size_up_cri = int(sys.argv[3])
	size_suff = size_suff + '_' + str(size_up_cri)
POS_col = 2-1
FILT_col = 7-1                                  # if this variant gets a 'PASS', it will be kept, without considerring the size
INFO_col = 8-1
ENDp = re.compile(r'(?<=\bEND=)(?P<end>\d+)(?=;|$)')
SVTYPEp = re.compile(r'(?<=\bSVTYPE=)(?P<svtype>[\w:\-_]+)(?=;|$)')
SVLENp =re.compile(r'(?<=\bSVLEN=)(?P<svlen>\-?\d+)(?=;|$)')

VCF=''
if vcf.endswith('.gz'):
	VCF=gzip.open(vcf,'rb')
else:
	VCF=file(vcf,'r')
annotations = []
header = ''
sv_filters = []

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
			pos = int(line.split('\t')[POS_col])
			end = int()
			svtype = ''
			svlen = int()
			try:
				svtype = SVTYPEp.search(line.split('\t')[INFO_col]).groupdict()['svtype']
				if svtype == 'BND':
				# translocation in Lumpy, ignore it
					continue
			except:
				sys.exit("Cannot find SVTYPE info in\n" + line )
			try:
				end = int(ENDp.search(line.split('\t')[INFO_col]).groupdict()['end'])
			except:
				try:
					svlen = int(SVLENp.search(line.split('\t')[INFO_col]).groupdict()['svlen'])
				except:
					sys.exit("1.Cannot find END or SVLEN info in\n" + line)
				end = pos if svlen == 1 else abs(svlen) + pos
			if svlen == 0:
				try:
					svlen = int(SVLENp.search(line.split('\t')[INFO_col]).groupdict()['svlen'])
				except:
					if end == 0:
						sys.exit("2.Cannot find END or SVLEN info in\n" + line)
					else:
						svlen = 1 if end == pos else end - pos
			size = 1 if end == pos else end - pos
			if(size != abs(svlen) and svlen != 1):
#				print size
#				print svlen
				print "You may check \n" + line
#			if size >= size_low_cri and size <= size_up_cri:
			if abs(svlen)>= size_low_cri and abs(svlen) <= size_up_cri:
				sv_filters.append(line)
VCF.close()

output = os.path.dirname(os.path.abspath(vcf)) + '/' + os.path.basename(vcf).strip('.gz').rpartition('.')[0] + '.' + size_suff + '.vcf'
OUT=file(output,'w')
OUT.writelines(annotations)
OUT.writelines(header)
OUT.writelines(sv_filters)
OUT.close()
