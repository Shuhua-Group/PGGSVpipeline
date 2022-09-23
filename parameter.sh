#software
PERL="/usr/bin/perl"
PYTHON="/usr/bin/python2.7"
SAMTOOLS="/share/apps/gene/samtools-1.14/bin/samtools"
BREAKDANCER="/home/wangyimin/anaconda3/envs/breakdancer/bin"
BREAKSEQ2="/home/wangyimin/anaconda3/envs/breakseq/bin"
BWA="/share/apps/gene/bwa-0.7.17/bwa"
CNVNATOR="/home/wangyimin/anaconda3/envs/cnvnator/bin"
LUMPY="/home/wangyimin/anaconda3/envs/lumpy/bin"
MANTA="/home/wangyimin/anaconda3/envs/manta/bin"
METASV="/home/wangyimin/anaconda3/envs/metasv/bin"

#parameters
reRunStop=5
BIN_SIZE=100
readlength=150
f1_size_l=1
f1_size_u=10000
f2_size_l=50
f2_size_u=2000000
RD_confSize=10000
CHR='all'
steps='all'
SOFT="all"
tools=(breakdancer breakseq2 cnvnator lumpy manta)
CHRlist="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"
HG19_CHRlist="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y"