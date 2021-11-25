#!/usr/bin/env python
## Used specifically to merge the DEL.nonOverlap and DUP.nonOverlap results, from MetaSV, after 'MetaSV_selfOverlap_v5.py' step.
## To be precise, I use 'CNVE' concept in this step. 
## 'CNV' means complex of DEL and DUP.

import os
import sys
import math
import re

ENDp=re.compile(r'(?<=\bEND=)(?P<end>\d+)(?=;|$)')
SVTYPEp = re.compile(r'(?<=\bSVTYPE=)(?P<svtype>[\w:\-_]+)(?=;|$)')
SOURCEp = re.compile(r'(?<=\bSOURCES=)(?P<source>[0-9A-Za-z\-,]+)(?=;|$)')

if len(sys.argv) != 3:
	sys.exit("Argument: MetaSV.DEL.nonOverlap.vcf MetaSV.DUP.nonOverlap.vcf")

DELvcf = sys.argv[1]
DUPvcf = sys.argv[2]
CHR_col = 1-1
POS_col = 2-1
FILT_col = 7-1
INFO_col = 8-1
verbose = 1000

Merge = []
annotations = []
header = ''

DEL = file(DELvcf,'r')
while True:
	line=DEL.readline().rstrip()
	if len(line) <= 0:
		break
	if line.startswith('##'):
		annotations.append(line + '\n')
	elif line.startswith('#'):
		try:
			if len(header) > 0:
				raise ExistedHeader
			else:
				header = line + '\n'
		except ExistedHeader:
			sys.exit('Multiple header in ' + MetaSV + '.')
	else:
		if(line[0] == 'X'):
			line = line.replace('X','23')
		if(line[0] == 'Y'):
			line = line.replace('Y','24')
		Merge.append(line.split())
DEL.close()
DUP = file(DUPvcf,'r')
while True:
	line=DUP.readline().rstrip()
	if len(line) <= 0:
		break
	if not line.startswith('#'):
		if(line[0] == 'X'):
			line = line.replace('X','23')
		if(line[0] == 'Y'):
			line = line.replace('Y','24')
		Merge.append(line.split())
DUP.close()

Merge.sort(key=lambda x:(int(x[0]),int(x[1])))

n=1
CNVE = []
while Merge:
	if verbose > 0:
		if n % verbose ==0:
			print str(n) + '...' 
	curLine = Merge.pop(0)
	n += 1
	curChr = curLine[CHR_col]
	curFilt = curLine[FILT_col]
	curPos = int(curLine[POS_col])
	curEnd = int()
	curType = ''
	curSource = ''
	try:
		curEnd = int(ENDp.search(curLine[INFO_col]).groupdict()['end'])
	except:
		sys.exit("1.No END info in line\n" + "\t".join(curLine) + "\n")
	try:
		curType = SVTYPEp.search(curLine[INFO_col]).groupdict()['svtype']
	except:
		sys.exit("1.No SVTYPE info in line\n" + "\t".join(curLine) + "\n")
	try:
		curSource = SOURCEp.search(curLine[INFO_col]).groupdict()['source']
	except:
		sys.exit("1.No SOURCE info in line\n" + "\t".join(curLine) + "\n")

	Overlaps = [[curLine,curEnd,curType]]			# link
	offset = curPos
	CNVstate = []
	posi2source = []
	posi2filter = []
	breakpoints = []
	for i in range((curPos-offset),(curEnd-offset+1)):
		if i != len(posi2source):
			print "1"
		if i != len(posi2filter):
			print "2"
		if i != len(CNVstate):
			print "3"
#		posi2source{i+offset} = curSource
#		posi2filter{i+offset} = curFilt
		posi2source.append(curSource)
		posi2filter.append(curFilt)
		CNVstate.append(curType)
#		try:
#			CNVstate[i] = 'CNV' if curType != CNVstate[i] else CNVstate[i]
#		except IndexError:
#			CNVstate[i] = curType
#		except:
#			sys.exit("maybe leading blank positions")
	tmp_end = Overlaps[-1][1]
	breakpoints=[curPos,curEnd]
	while Merge:
		tmp_chr = Overlaps[-1][0][CHR_col]
		tmp_pos = Overlaps[-1][0][POS_col]
		tmp_end = Overlaps[-1][1] if tmp_end <= Overlaps[-1][1] else tmp_end
		tmp_type = Overlaps[-1][2]
		nextLine = Merge[0]				# link
		nextChr = nextLine[CHR_col]
		nextFilt = nextLine[FILT_col]
		nextPos = int(nextLine[POS_col])
		nextEnd = int()
		nextType = ''
		nextSOurce = ''
		try:
			nextEnd = int(ENDp.search(nextLine[INFO_col]).groupdict()['end'])
		except:
			sys.exit("2.No END info in line\n" + "\t".join(nextLine) + "\n")
		try:
			nextType = SVTYPEp.search(nextLine[INFO_col]).groupdict()['svtype']
		except:
			sys.exit("2.No SVTYPE info in line\n" + "\t".join(nextLine) + "\n")
		try:
			nextSource = SOURCEp.search(nextLine[INFO_col]).groupdict()['source']
		except:
			sys.exit("2.No SOURCE info in line\n" + "\t".join(nextLine) + "\n")
		if nextChr == tmp_chr and nextPos <= tmp_end:
			Overlaps.append([Merge.pop(0),nextEnd,nextType])	# link
			n+=1
			breakpoints += [nextPos, nextEnd]
			for i in range((nextPos-offset),(nextEnd-offset+1)):
				if nextFilt == 'PASS':
					try:
						posi2filter[i] = 'PASS'
					except IndexError:
						if i != len(posi2filter):
							print "4"
						posi2filter.append('PASS')
					except:
						sys.exit("1.maybe leading blank positions")
				else:
					try:
						if posi2filter[i] != 'PASS':
							posi2filter[i] = nextFilt
					except IndexError:
						if i != len(posi2filter):
							print "5"
						posi2filter.append(nextFilt)
					except:
							sys.exit("2.maybe leading blank positions")
#				try:
#					posi2source{i+offset} += (',' + nextSource)
#				except KeyError:
#					posi2source{i+offset} = nextSource
#				except:
#					sys.exit("I don't know")
				try:
					posi2source[i] += (',' + nextSource)
				except IndexError:
					if i != len(posi2source):
						print "6"
					posi2source.append(nextSource)
				except:
					sys.exit("3.maybe leading blank positions")
				try:
					CNVstate[i] = 'CNV' if nextType != CNVstate[i] else CNVstate[i]
				except IndexError:
					if i != len(CNVstate):
						print "7"
					CNVstate.append(nextType)
				except:
					sys.exit("4.maybe leading blank positions")
		else:
			break
	if len(Overlaps) == 1:
		CNVE.append(Overlaps[0][0][:])
	else:
		BPs = sorted(set(breakpoints))
		for p in range(0,(len(BPs)-1)):
			start = BPs[p]
			end = BPs[p+1]
			SVLEN = end - start
			if SVLEN < 50:							# remove cnve less than 50bp
				continue
			cnveLine = [curChr,str(start),'.','N']
			midType = ''							# Ignore the CN state at start and end position, for some overlapping just at the breakpoints.
			states = set(CNVstate[(start+1-offset):(end-offset)])
			if 'CNV' in states:
				midType = 'CNV'
			elif ('DUP' in states) and ('DEL' in states):
				midType = 'CNV'
			elif len(states) == 1:
				midType = CNVstate[start+1-offset]
			else:
				sys.exit("Impossible!")
			ALT = '<' + midType + '>'
			FILTER = 'PASS' if 'PASS' in set(posi2filter[(start+1-offset):(end-offset)]) else 'LowQual'
			infos = []
			for i in list(set(posi2source[(start+1-offset):(end-offset)])):
				infos += i.split(',')
			NUM_SVTOOLS = len(set([meth.rpartition('-')[2] for meth in infos]))
			INFO = ";".join(["END="+str(end),"SVLEN="+str(SVLEN),"SVTYPE="+midType,"SVTOOL=MetaSV","SOURCES="+(','.join(infos)),"NUM_SVTOOLS="+str(NUM_SVTOOLS)])
			cnveLine += [ALT,'.',FILTER,INFO,'GT','1/1']
			CNVE.append(cnveLine)
	del(CNVstate,posi2source,posi2filter,breakpoints)

CNVEc = [("\t".join(cnve) + "\n") for cnve in CNVE]

CNVEvcf = os.path.dirname(os.path.abspath(DELvcf)) + '/' + os.path.basename(DELvcf).replace("DEL","DEL_DUP")
CNVEm = file(CNVEvcf,'w')
CNVEm.writelines(annotations)
CNVEm.writelines(header)
CNVEm.writelines(CNVEc)
CNVEm.close()
print "Finish!"

