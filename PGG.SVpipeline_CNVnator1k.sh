#!/bin/sh

##### Functions 

function usage0(){
	cat << INTRO
	[Usage]
	PGG.SVpipeline.sh --bam(f)ile indi.realigned.bam [[ --sample(i)d sampleid --(c)hr chromosome --work(d)ir working dir --(s)teps steps to run the pipeline --(m)etasvtools tools used to merge SVs by MetaSV --(r)eadlength readlength --(a)ssembly trun on assembly for MetaSV --(R)aw output raw results, instead of an embedded fitlering ] --dry specify running tools in dry-run model(, for tools already DONE) --(F)orce force to run MetaSV if any pre-metasv tool fails --(G)enotype perform genotyping after MetaSV step using CNVnator --gap_file genomic gaps to filter the CNV regions --gap_frac cutoff of gap filtering --pindel_vcf pindel vcf for MetaSV --breakdancer_naive breakdancer results for MetaSV --cnvnator_naive cnvnator results for MetaSV --breakseq2_vcf breakseq2 vcf for MetaSV --lumpy_vcf lumpy vcf for MetaSV ]
	-h for this short help
	--help for full help
INTRO
}
function usage1(){
	echo "-----------------------------------------------------------"
	echo "|in-house Pipeline for SV detection from NGS data, for PGG|"
	echo "-----------------------------------------------------------"
	cat << INTRO
	This script is used to call SVs/CNVs from realigned BAM file, to SV VCF file.
	The pipeline involves three components:
	1.) directly calling by different tools (i.e. different strategies), including: 
		CNVnator (RD)
		BreakDancer (RP)
		Lumpy (RP+SR)
		Pindel (SR)
		BreakSeq2 (JM)
	2.) merge (parts of) above results into single output, via:
		MetaSV
	3.) (optional) filtering for the results and output a final merged VCF, available for downstream analysis.

	[Usage]
		PGG.SVpipeline.sh -f|--bamfile individual.realigned.bam [options]
	[options]
		-f|--bamfile
			individual (realigned) BAM file [required]
		-i|--sampleid
			sample id; infer from bamfile, if not given
		-c|--chr
			chromosome to run; supported value: {1,2,...,22,23,24,X,Y,autosome,all}; 'all' as default
		-d|--workdir
			(top-level) working directory; sub-dir for each tool/sample will be made automatically; current dir will be uesd, if not given
		-s|--steps
			steps to run the pipeline; supported value: {cnvnator,lumpy,breakdancer,breakseq2,pindel,metasv,all}; multiple steps should be concatenated by ',' and order is not important; 'all' as default
		-m|--metasvtools
			(results from) tools used to merge SVs by MetaSV; supported value: {cnvnator,lumpy,breakdancer,breakseq2,pindel}; multiple tools should be concatenated by ',' and order is not important; all are used, if not supplied; required to supply 'metasv' to '--steps'
		-r|--readlength
			read length; 150 bp will be used, if not given
		-a|--assembly
			a switch to turn on the assembly for MetaSV (will be very slow); diabled as default
		-R|--Raw
			a switch to output very raw results; we have an embedded filtering process according to our prior knowledge to perform in the pipeline, as default.
		-F|--Force
			if any pre-metasv tool fails in one run (, even after re-run,) the MetaSV (if supplied in '--steps') step will not run as default, if this option is turn-on, the MetaSV will still run on the results of successfully run pre-meatsv tools
		-G|--Genotype
			perform an additional genotyping for SVs/CNVs from MetaSV merging results, using CNVnator; required to supply 'metasv' to '--steps'
		--gap_file
			regions to filter the calling CNVs (especially for DEL and DUP), in BED format, with first three columns as chr, start, and end. CNVs overlapping with these regions will be removed. Typically a gap file including centeromere and telomeres. (see --gap_frac)
		--gap_frac
			the cutoff of overlappings between CNV and gap regions (by --gap_file). larger than this value is considered as 'overlapped'. CNV size as denominator. (see --gap_file)

		(if you have existed results from pre-metasv tools and want to avoid running of such tools, use following options to supply the inputs)
		--pindel_vcf
			vcf results of pindel, supply to MetaSV
		--breakdancer_naive
			breakdancer results, supply to MetaSV
		--cnvnator_naive
			cnvnator results, supply to MetaSV
		--breakseq2_vcf
			vcf results of breakseq2, supply to MetaSV
		--lumpy_vcf
			vcf results of lumpy, supply to MetaSV

		(for special case, if you have finished (parts of) pre-metasv part of this pipeline, which means you have exactly the corresponding results in right format in terms of the filenames and the contents, you may use following '--dry' options to specify paticular step(s) to be run in 'dry-run' model (i.e. not really run, but treated as run, so that the existed results can be directly used for MetaSV step)

		--dry
			steps to run in dry-run model, especially that have been specified by '--steps' option; supported value: {cnvnator,lumpy,breakdancer,breakseq2,pindel}; multiple tools should be concatenated by ',' and order is not important;

			EXAMPLE:
				You first run the full pipeline with following command,
				>PGG.SVpipeline.sh -f Sample.realigned.bam
				but after some of the tools finished calling, say CNVnator, Pindel, and breakseq2, you have to kill the the process due to some unknown reasons. Then you fixed the problem and re-run the pipeline with following command,
				>PGG.SVpipeline.sh -f Sample.realigned.bam --steps all --metasvtools cnvnator,pindel,breakseq2,lumpy,breakdancer --dry cnvnator,pindel,breakseq2
				in this case, all steps will 'run' but cnvnator, pindel, and breakseq2 will run in 'dry-run' model, and then all the results (from five tools) are merged by MetaSV.
		
		-h 
			print short help
		--help
			print long help

	[Note]
	0.) this pipeline is especially useful to run from BAM files, throughout to get VCF files (in whole genome model). If you have existed pre-metasv results, we strongly suggest you to directly run MetaSV pipeline per se, otherwise we only provide single chromosome model for this situation! Or in a special situation, you may try '--dry' option.
	1.) although we provde 'chr' option, we only support 'whole genome' running for cnvnator, lumpy, and breakseq2.
	2.) breakdancer and pindel are run in divided chromosomes, even if you run in 'whole genome' model, considering the running time
	3.) all the runing parameters for involved tools are as default set, if you are an advanced user, you may modify the running command for any tool in the script, but you should be very clear about what you are doing.
	4.) there are some auto-check for failture running of some tools and will automatically re-run (till to success or server down) -- which is a double-edged sword
	5.) there is a also a stop marker (reRunStop) for 're-run' process, the default time is 5, which means if it continues fail, after 5 times of 're-run'-try, I will not try more, and you may manually check the error and then re-run the whole pipeline and you may consider using '--dry' option.

	[Example]
	1.) simple run
		PGG.SVpipeline.sh -f indi.realigned.bam 1>indi.SVpipeline.log 2>indi.SVpipeline.err 
		
	
	[Credit]
	Written by Fu Ruiqing, with contributions from Wu Zhendong, Lou Haiyi, and Zhang Xi.
	By courtesy of Tian Lei.

	[Contact]
	if you have any problem, please contact:
		Room 359
		furuiqing@picb.ac.cn

INTRO

}

function run_breakdancer(){
	export LD_LIBRARY_PATH=/picb/humpopg7/Frq_share/lib/:$LD_LIBRARY_PATH
	BAM=$1;	ID=$2; CHR=$3; WD=$4
	PATH=$(dirname $SAMTOOLS):$PATH $PERL $BREAKDANCER/perl/bam2cfg.pl $BAM > ${WD}/${ID}.breakdancer.cfg
	if [ "$CHR" = 'all' ] || [ "$CHR" = 'autosome' ]; then
		CHRs=$(echo {1..22}); [ "$CHR" = 'all' ] && CHRs=$(echo {1..22} X Y)
		for c in $CHRs
		do
			$BREAKDANCER/cpp/breakdancer_max -o $c ${WD}/${ID}.breakdancer.cfg > ${WD}/${ID}.chr${c}.SV.output &
		done
	elif [[ "$CHR" =~ ^[1-9]$ || "$CHR" =~ ^1[0-9]$ || "$CHR" =~ ^2[0-2]$ || "$CHR" =~ ^[XY]$ ]]; then
		$BREAKDANCER/cpp/breakdancer_max -o $CHR ${WD}/${ID}.breakdancer.cfg > ${WD}/${ID}.chr${CHR}.SV.output &
	else
		echo "BreakDancer: Unknown Chromosome!"; exit 1
	fi
	wait
	echo -n "$ID breakdancer DONE @ "; date
}

function run_lumpy(){
	BAM=$1;	ID=$2; CHR=$3; WD=$4
	if [ "$CHR" != 'all' ] && [ "$CHR" != 'autosome' ]; then
		touch ${WD}/.${ID}.lumpy.fail
		echo "We strongly suggest you to run Lumpy for Whole Genome!" ; exit 1
	fi
	$SAMTOOLS view -b -F 1294 $BAM > ${WD}/${ID}.discordants.unsorted.bam &
	$SAMTOOLS view -h $BAM | $LUMPY/scripts/extractSplitReads_BwaMem -i stdin | $SAMTOOLS view -Sb - > ${WD}/${ID}.splitters.unsorted.bam &
	wait
	$SAMTOOLS sort ${WD}/${ID}.discordants.unsorted.bam -o ${WD}/${ID}.discordants.bam &
	$SAMTOOLS sort ${WD}/${ID}.splitters.unsorted.bam -o ${WD}/${ID}.splitters.bam & # WYM edited syntax of samtools sort .bam -o .bam
	wait
	touch ${WD}/${ID}.lumpy.runlog
	grep "LUMPY Express done" ${WD}/${ID}.lumpy.runlog &>/dev/null && echo "Seems you have successfully ran Lumpy before. You may remove ${WD}/${ID}.lumpy.runlog file, if you are re-submitting for problems regarding Lumpy results."
	nrr=0
	until grep "LUMPY Express done" ${WD}/${ID}.lumpy.runlog &>/dev/null || [ "$nrr" -ge "$reRunStop" ]
	do
		if [ -s ${WD}/${ID}.lumpy.runlog ]; then
			echo "Error: re-run Lumpy for $ID!" >&2
			let nrr++
		fi
		$LUMPY/bin/lumpyexpress -B $BAM -S ${WD}/${ID}.splitters.bam -D ${WD}/${ID}.discordants.bam -o ${WD}/${ID}.lumpy.vcf > ${WD}/${ID}.lumpy.runlog &
		wait
	done
	# remove tmp files
	rm ${WD}/${ID}.discordants.unsorted.bam ${WD}/${ID}.splitters.unsorted.bam
	if grep "LUMPY Express done" ${WD}/${ID}.lumpy.runlog &>/dev/null; then
		echo -n "$ID lumpy DONE @ "; date
	else
		echo  -n "$ID lumpy FAILED @ "; date
		touch ${WD}/.${ID}.lumpy.fail
	fi
}

function run_pindel(){
	BAM=$1;	ID=$2; CHR=$3; WD=$4
	FAILf='FALSE'
	if [ "$CHR" = 'all' ] || [ "$CHR" = 'autosome' ]; then
		CHRs=$(echo {1..22}); [ "$CHR" = 'all' ] && CHRs=$(echo {1..22} X Y)
		echo "$BAM $INSize ${ID}" > ${WD}/${ID}.pindel.cfg
		for c in $CHRs
		do
			touch ${WD}/${ID}.chr${c}.runlog
			grep "PreviousChrName" ${WD}/${ID}.chr${c}.runlog &>/dev/null && echo "Seems you have successfully ran Pindel before. You may remove ${WD}/${ID}.chr${c}.runlog file, if you are re-submitting for problems regarding Pindel results."
			nrr=0
			until grep "PreviousChrName" ${WD}/${ID}.chr${c}.runlog &>/dev/null || [ "$nrr" -ge "$reRunStop" ]
			do
				if [ -s ${WD}/${ID}.chr${c}.runlog ]; then
					echo "Error: re-run Pindel for ${ID}.chr${c}!" >&2
					let nrr++
				fi
				$PINDEL/pindel -f $REF -T 2 -i ${WD}/${ID}.pindel.cfg -c $c -o ${WD}/${ID}.chr${c} > ${WD}/${ID}.chr${c}.runlog #louhaiyi edited -T 2
			done &
		done
		wait
		for c in $CHRs
		do
			if grep "PreviousChrName" ${WD}/${ID}.chr${c}.runlog &>/dev/null; then
				${PINDEL}/pindel2vcf -P ${WD}/${ID}.chr${c} -r $REF -R 1000GenomesPilot-NCBI37 -d 20090120 -v ${WD}/${ID}.pindel.chr${c}.vcf > ${WD}/${ID}.chr${c}.p2vlog &
			else
				FAILf='TRUE'
				echo -n "${ID}.chr${c} pindel FAILED @ ";date
			fi
		done
		wait
		[ "$FAILf" = 'FALSE' ] && echo -n "$ID Pindel DONE @ " && date
		[ "$FAILf" = 'TRUE' ] && touch ${WD}/.${ID}.pindel.fail
	elif [[ "$CHR" =~ ^[1-9]$ || "$CHR" =~ ^1[0-9]$ || "$CHR" =~ ^2[0-2]$ || "$CHR" =~ ^[XY]$ ]]; then
		echo "$BAM $INSize ${ID}.chr${CHR}" > ${WD}/${ID}.chr${CHR}.pindel.cfg
		touch ${WD}/${ID}.chr${CHR}.runlog
		grep "PreviousChrName" ${WD}/${ID}.chr${CHR}.runlog &>/dev/null && echo "Seems you have successfully ran Pindel before. You may remove ${WD}/${ID}.chr${CHR}.runlog file, if you are re-submitting for problems regarding Pindel results."
		nrr=0
		until grep "PreviousChrName" ${WD}/${ID}.chr${CHR}.runlog &>/dev/null || [ "$nrr" -ge "$reRunStop" ]
		do
			if [ -s ${WD}/${ID}.chr${CHR}.runlog ]; then
				echo "Error: re-run Pindel for ${ID}.chr${CHR}!" >&2
				let nrr++
			fi
			${PINDEL}/pindel -f $REF -i ${WD}/${ID}.chr${CHR}.pindel.cfg -c $CHR -o ${WD}/${ID}.chr${CHR} > ${WD}/${ID}.chr${CHR}.runlog
		done
		if grep "PreviousChrName" ${WD}/${ID}.chr${CHR}.runlog &>/dev/null; then
			${PINDEL}/pindel2vcf -P ${WD}/${ID}.chr${CHR} -r $REF -R 1000GenomesPilot-NCBI37 -d 20090120 -v ${WD}/${ID}.pindel.chr${CHR}.vcf > ${WD}/${ID}.chr${CHR}.p2vlog
			echo -n "$ID Pindel DONE @ "; date
		else
			echo -n "${ID}.chr${CHR} pindel FAILED @ "; date
			touch ${WD}/.${ID}.pindel.fail
		fi
	else
		echo "Pindel: Unknown Chromosome!"; exit 1
	fi
#	echo -n "$ID Pindel DONE @ "; date
}

function run_cnvnator(){
	source /picb/humpopg7/Frq_share/ROOT/bin/thisroot.sh
	BAM=$1;	ID=$2; CHR=$3; WD=$4
	if [ "$CHR" != 'all' ] && [ "$CHR" != 'autosome' ]; then
		touch ${WD}/.${ID}.cnvnator.fail
		echo "We strongly suggest you to run cnvnator for Whole Genome!"; exit 1
	fi
	touch ${WD}/${ID}.cnvnator.runlog
	$CNVNATOR -root ${WD}/${ID}.root -unique -genome GRCh37 -tree $BAM > ${WD}/${ID}.cnvnator.runlog
	$CNVNATOR -root ${WD}/${ID}.root -genome GRCh37 -his $BIN_SIZE -d $REFc > ${WD}/${ID}.cnvnator.runlog
	$CNVNATOR -root ${WD}/${ID}.root -stat $BIN_SIZE > ${WD}/${ID}.cnvnator.runlog
	$CNVNATOR -root ${WD}/${ID}.root -partition $BIN_SIZE > ${WD}/${ID}.cnvnator.runlog
	$CNVNATOR -root ${WD}/${ID}.root -call $BIN_SIZE > ${WD}/${ID}.rawcnv
	
	# for the format of running meatsv
	sed -i -e 's/\bchr//' ${WD}/${ID}.rawcnv
	echo -n "$ID cnvnator DONE @ "; date
}

function run_breakseq2(){
	BAM=$1;	ID=$2; CHR=$3; WD=$4
	if [ "$CHR" != 'all' ] && [ "$CHR" != 'autosome' ]; then
		touch ${WD}/.${ID}.breakseq2.fail
		echo "We strongly suggest you to run breakseq2 for Whole Genome!"; exit 1
	fi
	touch ${WD}/${ID}.breakseq2.runlog
	grep "breakseq2_workflow   Done" ${WD}/${ID}.breakseq2.runlog &>/dev/null && echo "Seems you have successfully ran breakseq2 before. You may remove ${WD}/${ID}.breakseq2.runlog file, if you are re-submitting for problems regarding breakseq2 results."
	nrr=0
	until grep "breakseq2_workflow   Done" ${WD}/${ID}.breakseq2.runlog &>/dev/null || [ "$nrr" -ge "$reRunStop" ]
	do
		if [ -s ${WD}/${ID}.breakseq2.runlog ]; then
			echo "Error: re-run breakseq2 for $ID!" >&2
			let nrr++
		fi
		$BREAKSEQ2 --nthreads 1 --bams $BAM --work ${WD} --reference $REF --sample $ID --bplib $BS2_bplib --bwa $BWA --samtools $SAMTOOLS --min_span 10 --window $BIN_SIZE --min_overlap 10 2>${WD}/${ID}.breakseq2.runlog
	done
	if grep "breakseq2_workflow   Done" ${WD}/${ID}.breakseq2.runlog &>/dev/null; then
		echo -n "$ID breakseq2 DONE @ "; date
	else
		echo -n "$ID breakseq2 FAILED @ "; date
		touch ${WD}/.${ID}.breakseq2.fail
	fi
}


##### set environment and TOOL/ANNOTATION path
which perl &>/dev/null || exit 1
which python &>/dev/null || exit 1
export BIN=$(dirname $0)	# script directory
export BIN_SIZE=100		# bin/win size for read-depth methods
export INSize=400		# insert size
export INSd=50			# insert sd
export REFc=/picb/humpopg-bigdata5/wangyimin/resource/Human_Genome_hg19/GRCh37_fa_files
export REF=/picb/humpopg-bigdata5/wangyimin/resource/Human_Genome_hg19/human_g1k_v37.fasta
export BS2_bplib=/picb/humpopg-bigdata/furuiqing/Resource/Breakpoints/breakseq2_bplib_20150129.fna
export PERL=/picb/humpopg7/Frq_share/Perl/bin/perl
export PYTHON=/picb/humpopg7/Frq_share/bin/python
export BWA=/picb/pggtools/archives/bwa-0.7.17/bwa
export SAMTOOLS=/home/wangyimin/samtools-1.3.1/samtools
export PINDEL=/picb/pggtools/archives/pindel-0.2.5b8
export BREAKSEQ2=/picb/humpopg7/Frq_share/python/bin/run_breakseq2.py
export CNVNATOR=/picb/humpopg7/Frq_share/CNVnator/CNVnator_v0.3/src/cnvnator
export BREAKDANCER=/picb/humpopg-bigdata3/wangyimin/Softs/breakdancer-1.1_2011_02_21
export LUMPY=/picb/pggtools/archives/lumpy-sv-0.2.13
export METASV=/picb/humpopg7/Frq_share/MetaSV/metasv-0.5/MSpy/bin/run_metasv.py
export SPADES=/picb/humpopg7/Frq_share/SPAdes/SPAdes-3.6.1-Linux/bin/spades.py
export AGE=/picb/humpopg7/Frq_share/AGE/AGE-master/age_align
###################################

##### embedded filter parameters
# first for non-cnvnator results
export f1_size_l=1
export f1_size_u=10000
# second for metasv results
export f2_size_l=50
export f2_size_u=2000000
export RD_confSize=1000
###################################


##### set variables
BAM=''
ID=''
CHR=''
WD=$PWD
steps='all'
geno='FALSE'
force='FALSE'
gap_frac=0.0
gap_file="NULL"
TOOLs=(breakdancer lumpy pindel cnvnator breakseq2 metasv)
export reRunStop=5
###################################

# use getopt to get all the arguments
GETOPT_ARGS=`getopt -n 'argument error' -o f:i:c:d:s:m:r:aRGFh -al bamfile:,sampleid:,chr:,workdir:,steps:,metasvtools:,dry:,Genotype,Force,gap_file:,gap_frac:,pindel_vcf:,breakdancer_naive:,cnvnator_naive:,breakseq2_vcf:,lumpy_vcf:,readlength:,assembly,Raw,help  -- "$@"` ||  eval "usage0; exit 1"
eval set -- "$GETOPT_ARGS"

# assign arguments
echo "###########Paramter setting#########"
while [ -n "$1" ] 
do
	case "$1" in 
		-f|--bamfile) BAM=$2; echo "bamfile: $2"; shift 2;;
		-i|--sampleid) ID=$2; echo "sampleid: $2"; shift 2;;
		-c|--chr) CHR=$2; echo "chromosome: $2"; shift 2;;
		-d|--workdir) WD=$2; echo "workdir: $2"; shift 2;;
		--workDir) WD=$2; echo "workdir: $2"; shift 2;;
		-s|--steps) steps=$2; echo "steps: $2"; shift 2;;
		-m|--metasvtools) MSts=$2; echo "tools for MetaSV: $2"; shift 2;;
		--dry) dry=$2; echo "dry run on: $2"; shift 2;;
		-r|--readlength) readlength=$2; echo "read length: $2"; shift 2;;
		-a|--assembly) MS_assembly="TRUE"; echo "will run assembly of MetaSV"; shift;;
		-R|--Raw) Raw="TRUE"; echo "will NOT perform the embedded filter process"; shift;;
		-G|--Genotype) geno='TRUE'; echo "will perform genotyping on CNVRs using CNVnator"; shift;;
		-F|--Force) force='TRUE'; echo "will run MetaSV, if you insist"; shift;;
		--gap_file) gap_file=$2; echo "gap file: $2"; shift 2;;
		--gap_frac) gap_frac=$2; echo "gap overlapping cutoff to filt CNV regions (fraction as CNV): >$2"; shift 2;;
		--pindel_vcf) pindelvcf=$2; echo "existed pindel vcf: $2"; shift 2;;
		--breakdancer_naive) breakdancerout=$2; echo "existed breakdancer output: $2"; shift 2;;
		--cnvnator_naive) cnvnatorout=$2; echo "existed cnvnator output: $2"; shift 2;;
		--breakseq2_vcf) breakseq2vcf=$2; echo "existed breakseq2 vcf: $2"; shift 2;;
		--lumpy_vcf) lumpyvcf=$2; echo "existed lumpy vcf: $2"; shift 2;;
		-h) usage0; exit 1;;
		--help) usage1; exit 1;;
		--) shift; break;;
		*) echo  "Unknown argument $1" ;usage1; exit 1;;
	esac
done

if [ -z "$BAM" ] || ! [ -e "$BAM" ];then
	echo "No available BAM file $BAM"
	echo "####################################"
	usage0; exit 1
fi

if [ -z "$ID" ]; then
	ID=$(basename $BAM .bam); ID=$(basename $ID .realigned)
	echo "You are not giving me a sampleid, I guess it is $ID"
fi

if [ -z "$CHR" ]; then
	echo "No specified chromosome, will run on whole genome!"
	CHR='all'
fi
CHR=$(echo $CHR |tr '[:upper:]' '[:lower:]')
CHR=${CHR#chr}
[ "$CHR" = 'x' -o "$CHR" = '23' ] && CHR='X'
[ "$CHR" = 'y' -o "$CHR" = '24' ] && CHR='Y'

if [ -z "$readlength" ]; then
	readlength=150
	echo "You are not giving me the readlength, I guess it is ${readlength}bp"
fi

if [ -z "$MS_assembly" ]; then
	MS_assembly='--disable_assembly'
else
	MS_assembly="--isize_mean $INSize --isize_sd $INSd --bam ${BAM} --spades $SPADES --age $AGE --boost_sc"
fi

if [ -z "$Raw" ]; then
	Raw="FALSE"
fi

# get steps 
steps=$(echo $steps |tr '[:upper:]' '[:lower:]')
STEPS=(${steps//,/ })
run=()
echo "You are going to run following tools/pipeline: "
for step in ${STEPS[@]}
do
	if [ "$step" = "all" ]; then
		run=${TOOLs[@]}
		echo -n "Ok, will run the entire pipeline @ "; date
		break
	elif echo "${TOOLs[*]}" |grep -w "$step" &>/dev/null; then
		echo "$step"
		run+=($step)
	else
		echo "$step --> OMG!! It seems that we currently don't support tool $step, please check if it is a typo!"
		echo "####################################"
		usage0; exit 1
	fi
done

# get dry-runs
dry=$(echo $dry |tr '[:upper:]' '[:lower:]')
dry=(${dry//,/ })
norun=()
for step in ${dry[@]}
do
	if echo "${run[*]}" |grep -w "$step" &>/dev/null; then
		echo "will do dry-run on $step"
		norun+=($step)
	else
		echo "$step is not going to be running, ignore"
	fi
done

# get metasvtools
if [ -z "$MSts" ]; then
	if echo "${run[*]}" |grep -w "metasv" &>/dev/null; then
		MSts=(${run[*]//metasv/})
		if [ ${#MSts[@]} -eq 0 ]; then
			MSts=${TOOLs[@]}; MSts=(${MSts[*]//metasv/})
		fi
	fi
else
	if ! echo "${run[*]}" |grep -w "metasv" &>/dev/null; then
		echo "You are not asking me to run MetaSV, but why give me 'metasvtools'"
		echo "####################################"
		usage0; exit 1
	fi
	MSts=$(echo $MSts |tr '[:upper:]' '[:lower:]')
	MSts=(${MSts//,/ })
fi


if [ -n "$MSts" ]; then

	echo "You are going to merge results from following tools (by MetaSV): "
	echo ${MSts[*]} |tr -s ' ' ','

	echo "${MSts[*]}" | grep -w "lumpy" &>/dev/null && (echo "${run[*]}" |grep -w "lumpy" &>/dev/null && [ -n "$lumpyvcf" ]) && echo "You don't need to provide Lumpy vcf for MetaSv if you are going to run Lumpy at the same time!" && exit 1
	echo "${MSts[*]}" | grep -w "pindel" &>/dev/null && (echo "${run[*]}" |grep -w "pindel" &>/dev/null && [ -n "$pindelvcf" ]) && echo "You don't need to provide Pindel vcf for MetaSv if you are going to run Pindel at the same time!" && exit 1
	echo "${MSts[*]}" | grep -w "cnvnator" &>/dev/null && (echo "${run[*]}" |grep -w "cnvnator" &>/dev/null && [ -n "$cnvnatorout" ]) && echo "You don't need to provide CNVnator output  for MetaSv if you are going to run CNVnator at the same time!" && exit 1
	echo "${MSts[*]}" | grep -w "breakdancer" &>/dev/null && (echo "${run[*]}" |grep -w "breakdancer" &>/dev/null && [ -n "$breakdancerout" ]) && echo "You don't need to provide BreakDancer output for MetaSv if you are going to run BreakDancer at the same time!" && exit 1
	echo "${MSts[*]}" | grep -w "breakseq2" &>/dev/null && (echo "${run[*]}" |grep -w "breakseq2" &>/dev/null && [ -n "$breakseq2vcf" ]) && echo "You don't need to provide BreakSeq2 vcf for MetaSv if you are going to run BreakSeq2 at the same time!" && exit 1
	echo "${MSts[*]}" | grep -w "lumpy" &>/dev/null && (! echo "${run[*]}" |grep -w "lumpy" &>/dev/null && ! [ -n "$lumpyvcf" ]) && echo "You are neither running Lumpy nor providing me a Lumpy vcf to run MetaSV!" && exit 1
	echo "${MSts[*]}" | grep -w "pindel" &>/dev/null && (! echo "${run[*]}" |grep -w "pindel" &>/dev/null && ! [ -n "$pindelvcf" ]) && echo "You are neither running Pindel nor providing me a Pindel vcf to run MetaSV!" && exit 1
	echo "${MSts[*]}" | grep -w "cnvnator" &>/dev/null && (! echo "${run[*]}" |grep -w "cnvnator" &>/dev/null && ! [ -n "$cnvnatorout" ]) && echo "You are neither running CNVnator nor providing me a CNVnator output to run MetaSV!" && exit 1
	echo "${MSts[*]}" | grep -w "breakdancer" &>/dev/null && (! echo "${run[*]}" |grep -w "breakdancer" &>/dev/null && ! [ -n "$breakdancerout" ]) && echo "You are neither running BreakDancer nor providing me a BreakDancer output to run MetaSV!" && exit 1
	echo "${MSts[*]}" | grep -w "breakseq2" &>/dev/null && (! echo "${run[*]}" |grep -w "breakseq2" &>/dev/null && ! [ -n "$breakseq2vcf" ]) && echo "You are neither running BreakSeq2 nor providing me a BreakSeq2 vcf to run MetaSV!" && exit 1
fi

[ "$geno" = "TRUE" ] && ! echo "${run[*]}" |grep -w "metasv" &>/dev/null && echo "additional genotyping is provided for after running of MetaSV, but you are not going to run MetaSV!" && exit 1

echo "####################################"	# parameter setting

## Start to run
run_MetaSV='FALSE'
premetasv_flag=(TRUE TRUE TRUE TRUE TRUE)	# run success for pre-metasv tools, (breakdancer, lumpy, pindel, cnvnator, breakseq2)
for tool in ${run[@]}
do
	if [ "$tool" = 'breakdancer' ]; then
		[ ! -d "${WD}/breakdancer/${ID}" ] && mkdir -p -m 775 "${WD}/breakdancer/${ID}"
		if echo "${norun[*]}" |grep -w "breakdancer" &>/dev/null; then
			echo -n "$ID breakdancer (dry-run) DONE @ "; date
		else
			run_breakdancer $BAM $ID $CHR ${WD}/breakdancer/${ID} 2>${WD}/${ID}.breakdancer.err &
		fi
	elif [ "$tool" = 'lumpy' ]; then
		[ ! -d "${WD}/lumpy/${ID}" ] && mkdir -p -m 775 "${WD}/lumpy/${ID}" 
		if echo "${norun[*]}" |grep -w "lumpy" &>/dev/null; then
			echo -n "$ID lumpy (dry-run) DONE @ "; date
		else
			run_lumpy $BAM $ID $CHR ${WD}/lumpy/${ID} 2>${WD}/${ID}.lumpy.err &
		fi
	elif [ "$tool" = 'pindel' ]; then
		[ ! -d "${WD}/pindel/${ID}" ] && mkdir -p -m 775 "${WD}/pindel/${ID}" 
		if echo "${norun[*]}" |grep -w "pindel" &>/dev/null; then
			echo -n "$ID pindel (dry-run) DONE @ "; date
		else
			run_pindel $BAM $ID $CHR ${WD}/pindel/${ID} 2>${WD}/${ID}.pindel.err &
		fi
	elif [ "$tool" = 'cnvnator' ]; then
		[ ! -d "${WD}/cnvnator/${ID}" ] && mkdir -p -m 775 "${WD}/cnvnator/${ID}" 
		if echo "${norun[*]}" |grep -w "cnvnator" &>/dev/null; then
			echo -n "$ID cnvnator (dry-run) DONE @ "; date
		else
			run_cnvnator $BAM $ID $CHR ${WD}/cnvnator/${ID} 2>${WD}/${ID}.cnvnator.err &
		fi
	elif [ "$tool" = 'breakseq2' ]; then
		[ ! -d "${WD}/breakseq2/${ID}" ] && mkdir -p -m 775 "${WD}/breakseq2/${ID}"
		if echo "${norun[*]}" |grep -w "breakseq2" &>/dev/null; then
			echo -n "$ID breakseq2 (dry-run) DONE @ "; date
		else
			run_breakseq2 $BAM $ID $CHR ${WD}/breakseq2/${ID} 2>${WD}/${ID}.breakseq2.err &
		fi
	elif [ "$tool" = 'metasv' ]; then
		[ ! -d "${WD}/metasv/${ID}" ] && mkdir -p -m 775 "${WD}/metasv/${ID}"
		if echo "${norun[*]}" |grep -w "metasv" &>/dev/null; then
			echo "We do NOT suggest run MetaSV in dry-run model. Anyway, I will continue this time (by skip MetaSV step), but be cautious next time!"
		else
			run_MetaSV='TRUE'
		fi
	else
		echo "I don't know the tool $tool! Quit now..."; exit 1
	fi
done
wait

## Check failtures
for tool in ${MSts[@]}
do
	[ "$tool" = "breakdancer" ] && [ -e "${WD}/breakdancer/${ID}/.${ID}.breakdancer.fail" ] && premetasv_flag[0]='FALSE' && rm "${WD}/breakdancer/${ID}/.${ID}.breakdancer.fail"
	[ "$tool" = "lumpy" ] && [ -e "${WD}/lumpy/${ID}/.${ID}.lumpy.fail" ] && premetasv_flag[1]='FALSE' && rm "${WD}/lumpy/${ID}/.${ID}.lumpy.fail"
	[ "$tool" = "pindel" ] && [ -e "${WD}/pindel/${ID}/.${ID}.pindel.fail" ] && premetasv_flag[2]='FALSE' && rm "${WD}/pindel/${ID}/.${ID}.pindel.fail"
	[ "$tool" = "cnvnator" ] && [ -e "${WD}/cnvnator/${ID}/.${ID}.cnvnator.fail" ] && premetasv_flag[3]='FALSE' && rm "${WD}/cnvnator/${ID}/.${ID}.cnvnator.fail"
	[ "$tool" = "breakseq2" ] && [ -e "${WD}/breakseq2/${ID}/.${ID}.breakseq2.fail" ] && premetasv_flag[4]='FALSE' && rm "${WD}/breakseq2/${ID}/.${ID}.breakseq2.fail"
done

## run embedded filter process and merge by MetaSV
filterName=''
MetaSVf='TRUE'
if [ "$run_MetaSV" = 'FALSE' ]; then
	echo -n "Finish all callings (, but MetaSV) @ "; date
elif [ "$force" = 'FALSE' ] && echo "${premetasv_flag[*]}" |grep -w "FALSE" &>/dev/null; then
	echo "Error: Stop MetaSV merging, due to some pre-metasv tools failing to run successfully, check above logs for help!"
	echo -n "Abandon pipeline @ "; date
	MetaSVf='FALSE'
else
	echo "${premetasv_flag[*]}" |grep -w "FALSE" &>/dev/null && echo "Some pre-metasv tools does NOT successfully run, continue merging without them..."
	echo -n "Start running MetaSV @ "; date

	if [ "$CHR" = 'all' ] || [ "$CHR" = 'autosome' ]; then
		CHRs=$(echo {1..22}); [ "$CHR" = 'all' ] && CHRs=$(echo {1..22} X Y)
		[ ! -d "${WD}/metasv/${ID}" ] && mkdir -p -m 775 "${WD}/metasv/${ID}" 
		MetaSV_arg_ini="--reference ${REF} --sample ${ID} --num_threads 1 --enable_per_tool_output --keep_standard_contigs --mean_read_length $readlength $MS_assembly"
		for c in $CHRs
		do
			MetaSV_arg="$MetaSV_arg_ini --outdir ${WD}/metasv/${ID}/chr${c} --workdir ${WD}/metasv/${ID}/chr${c}/tmp_work --chromosomes $c"
			if  echo "${MSts[*]}" | grep -w "lumpy" &>/dev/null && [ "${premetasv_flag[1]}" = 'TRUE' ]; then
				lumpyinput=${lumpyvcf:-${WD}/lumpy/${ID}/${ID}.lumpy.vcf}
				if [ "$Raw" = 'FALSE' ]; then
					filterName=".size_${f1_size_l}_${f1_size_u}"
#					$PYTHON $BIN/VCFfilterSize.py ${lumpyinput} $f1_size_l $f1_size_u > ${WD}/lumpy/${ID}/${ID}.lumpy.runfilter.log
					[ ! -d "${WD}/lumpy" ] && mkdir -p -m 775 "${WD}/lumpy"
					[ -d "${WD}/lumpy/${ID}/" ] && $PYTHON $BIN/VCFfilterSize.py ${lumpyinput} $f1_size_l $f1_size_u > ${WD}/lumpy/${ID}/${ID}.lumpy.runfilter.log || $PYTHON $BIN/VCFfilterSize.py ${lumpyinput} $f1_size_l $f1_size_u > ${WD}/lumpy/${ID}.lumpy.runfilter.log
					t=$(basename $lumpyinput .gz); lumpyinput="$(dirname $lumpyinput)/${t%.*}${filterName}.vcf"
				fi
				MetaSV_arg="$MetaSV_arg --lumpy_vcf $lumpyinput"
			fi
			if  echo "${MSts[*]}" | grep -w "cnvnator" &>/dev/null && [ "${premetasv_flag[3]}" = 'TRUE' ]; then
				MetaSV_arg="$MetaSV_arg --cnvnator_native ${cnvnatorout:-${WD}/cnvnator/${ID}/${ID}.rawcnv}"
			fi
			if  echo "${MSts[*]}" | grep -w "breakdancer" &>/dev/null && [ "${premetasv_flag[0]}" = 'TRUE' ]; then
				breakdancerinput=${breakdancerout:-${WD}/breakdancer/${ID}/${ID}.chr${c}.SV.output}
				if [ "$Raw" = 'FALSE' ]; then
					filterName=".size_${f1_size_l}_${f1_size_u}"
#					$BIN/BreakDancerSizeFilter.sh ${breakdancerinput} $f1_size_l $f1_size_u > ${WD}/breakdancer/${ID}/${ID}.chr${c}.runfilter.log
					[ ! -d "${WD}/breakdancer" ] && mkdir -p -m 775 "${WD}/breakdancer"
					[ -d "${WD}/breakdancer/${ID}" ] && $BIN/BreakDancerSizeFilter.sh ${breakdancerinput} $f1_size_l $f1_size_u > ${WD}/breakdancer/${ID}/${ID}.chr${c}.runfilter.log || $BIN/BreakDancerSizeFilter.sh ${breakdancerinput} $f1_size_l $f1_size_u > ${WD}/breakdancer/${ID}.chr${c}.runfilter.log
					t=$(basename $breakdancerinput .SV.output); breakdancerinput="$(dirname $breakdancerinput)/${t}${filterName}.SV.output"
				fi
				MetaSV_arg="$MetaSV_arg --breakdancer_native $breakdancerinput"
			fi
			if  echo "${MSts[*]}" | grep -w "breakseq2" &>/dev/null && [ "${premetasv_flag[4]}" = 'TRUE' ]; then
				breakseq2input=${breakseq2vcf:-${WD}/breakseq2/${ID}/breakseq.vcf.gz}
				if [ "$Raw" = 'FALSE' ]; then
					filterName=".size_${f1_size_l}_${f1_size_u}"
#					$PYTHON $BIN/VCFfilterSize.py ${breakseq2input} $f1_size_l $f1_size_u > ${WD}/breakseq2/${ID}/${ID}.runfilter.log
					[ ! -d "${WD}/breakseq2" ] && mkdir -p -m 775 "${WD}/breakseq2"
					[ -d "${WD}/breakseq2/${ID}" ] && $PYTHON $BIN/VCFfilterSize.py ${breakseq2input} $f1_size_l $f1_size_u > ${WD}/breakseq2/${ID}/${ID}.runfilter.log || $PYTHON $BIN/VCFfilterSize.py ${breakseq2input} $f1_size_l $f1_size_u > ${WD}/breakseq2/${ID}.runfilter.log
					t=$(basename $breakseq2input .gz); breakseq2input="$(dirname $breakseq2input)/${t%.*}${filterName}.vcf"
				fi
				MetaSV_arg="$MetaSV_arg --breakseq_vcf $breakseq2input"
			fi
			if  echo "${MSts[*]}" | grep -w "pindel" &>/dev/null && [ "${premetasv_flag[2]}" = 'TRUE' ]; then
				pindelinput=${pindelvcf:-${WD}/pindel/${ID}/${ID}.pindel.chr${c}.vcf}
				if [ "$Raw" = 'FALSE' ]; then
					filterName=".size_${f1_size_l}_${f1_size_u}"
#					$PYTHON $BIN/VCFfilterSize.py $pindelinput $f1_size_l $f1_size_u > ${WD}/pindel/${ID}/${ID}.pindel.chr${c}.runfilter.log
					[ ! -d "${WD}/pindel" ] && mkdir -p -m 775 "${WD}/pindel"
					[ -d "${WD}/pindel/${ID}" ] && $PYTHON $BIN/VCFfilterSize.py $pindelinput $f1_size_l $f1_size_u > ${WD}/pindel/${ID}/${ID}.pindel.chr${c}.runfilter.log || $PYTHON $BIN/VCFfilterSize.py $pindelinput $f1_size_l $f1_size_u > ${WD}/pindel/${ID}.pindel.chr${c}.runfilter.log
					t=$(basename $pindelinput .gz); pindelinput="$(dirname $pindelinput)/${t%.*}${filterName}.vcf"
				fi
				MetaSV_arg="$MetaSV_arg --pindel_vcf $pindelinput"
			fi
			touch ${WD}/metasv/${ID}/${ID}.metasv.chr${c}.err
			grep "metasv\.main\s*All Done" ${WD}/metasv/${ID}/${ID}.metasv.chr${c}.err &>/dev/null && echo "Seems you have successfully ran metasv before. You may remove ${WD}/${ID}/${ID}.metasv.chr${c}.err file, if you are re-submitting for problems regarding metasv results."
			nrr=0
			until grep "metasv\.main\s*All Done"  ${WD}/metasv/${ID}/${ID}.metasv.chr${c}.err &>/dev/null || [ "$nrr" -ge "$reRunStop" ]
			do
				if [ -s ${WD}/metasv/${ID}/${ID}.metasv.chr${c}.err ]; then
					echo "Error: re-run MetaSV for $ID chr${c}!"
					let nrr++
				fi
				$METASV $MetaSV_arg 2>${WD}/metasv/${ID}/${ID}.metasv.chr${c}.err
			done &
		done
		wait
		FAILf='FALSE'
		for c in $CHRs
		do
			! grep "metasv\.main\s*All Done"  ${WD}/metasv/${ID}/${ID}.metasv.chr${c}.err &>/dev/null && FAILf='TRUE' && echo "MetaSV failed for chromosome $c"
		done
		if [ "$FAILf" = 'FALSE' ]; then
			$BIN/MergeIndiChr.pl ${WD}/metasv/${ID} $ID > ${WD}/metasv/${ID}/${ID}.merge.log
		else
			echo "MetaSV failed for some chromosomes, will NOT merge the entire genome!"
			MetaSVf='FALSE'
		fi
	elif [[ "$CHR" =~ ^[1-9]$ || "$CHR" =~ ^1[0-9]$ || "$CHR" =~ ^2[0-2]$ || "$CHR" =~ ^[XY]$ ]]; then
		[ ! -d "${WD}/metasv/${ID}/chr${CHR}" ] && mkdir -p -m 775 "${WD}/metasv/${ID}/chr${CHR}" 
		MetaSV_arg="--reference ${REF} --outdir ${WD}/metasv/${ID}/chr${CHR} --sample ${ID} --num_threads 1 --workdir ${WD}/metasv/${ID}/chr${CHR}/tmp_work --enable_per_tool_output --keep_standard_contigs --chromosomes $CHR --mean_read_length $readlength $MS_assembly"
		if  echo "${MSts[*]}" | grep -w "lumpy" &>/dev/null && [ "${premetasv_flag[1]}" = 'TRUE' ]; then
			lumpyinput=${lumpyvcf:-${WD}/lumpy/${ID}/${ID}.lumpy.vcf}
			if [ "$Raw" = 'FALSE' ]; then
				filterName=".size_${f1_size_l}_${f1_size_u}"
#				$PYTHON $BIN/VCFfilterSize.py $lumpyinput $f1_size_l $f1_size_u > ${WD}/lumpy/${ID}/${ID}.lumpy.runfilter.log
				[ ! -d "${WD}/lumpy" ] && mkdir -p -m 775 "${WD}/lumpy"
				[ -d "${WD}/lumpy/${ID}" ] && $PYTHON $BIN/VCFfilterSize.py $lumpyinput $f1_size_l $f1_size_u > ${WD}/lumpy/${ID}/${ID}.lumpy.runfilter.log || $PYTHON $BIN/VCFfilterSize.py $lumpyinput $f1_size_l $f1_size_u > ${WD}/lumpy/${ID}.lumpy.runfilter.log
				t=$(basename $lumpyinput .gz); lumpyinput="$(dirname $lumpyinput)/${t%.*}${filterName}.vcf"
			fi
			MetaSV_arg="$MetaSV_arg --lumpy_vcf $lumpyinput"
		fi
		if  echo "${MSts[*]}" | grep -w "cnvnator" &>/dev/null && [ "${premetasv_flag[3]}" = 'TRUE' ]; then
			MetaSV_arg="$MetaSV_arg --cnvnator_native ${cnvnatorout:-${WD}/cnvnator/${ID}/${ID}.rawcnv}"
		fi
		if  echo "${MSts[*]}" | grep -w "breakdancer" &>/dev/null && [ "${premetasv_flag[0]}" = 'TRUE' ]; then
			breakdancerinput=${breakdancerout:-${WD}/breakdancer/${ID}/${ID}.chr${CHR}.SV.output}
			if [ "$Raw" = 'FALSE' ]; then
				filterName=".size_${f1_size_l}_${f1_size_u}"
#				$BIN/BreakDancerSizeFilter.sh $breakdancerinput $f1_size_l $f1_size_u > ${WD}/breakdancer/${ID}/${ID}.chr${CHR}.runfilter.log
				[ ! -d "${WD}/breakdancer" ] && mkdir -p -m 775 "${WD}/breakdancer"
				[ -d "${WD}/breakdancer/${ID}" ] && $BIN/BreakDancerSizeFilter.sh $breakdancerinput $f1_size_l $f1_size_u > ${WD}/breakdancer/${ID}/${ID}.chr${CHR}.runfilter.log || $BIN/BreakDancerSizeFilter.sh $breakdancerinput $f1_size_l $f1_size_u > ${WD}/breakdancer/${ID}.chr${CHR}.runfilter.log
				t=$(basename $breakdancerinput .SV.output); breakdancerinput="$(dirname $breakdancerinput)/${t}${filterName}.SV.output"
			fi
			MetaSV_arg="$MetaSV_arg --breakdancer_native $breakdancerinput"
		fi
		if  echo "${MSts[*]}" | grep -w "breakseq2" &>/dev/null && [ "${premetasv_flag[4]}" = 'TRUE' ]; then
			breakseq2input=${breakseq2vcf:-${WD}/breakseq2/${ID}/breakseq.vcf.gz}
			if [ "$Raw" = 'FALSE' ]; then
				filterName=".size_${f1_size_l}_${f1_size_u}"
#				$PYTHON $BIN/VCFfilterSize.py $breakseq2input $f1_size_l $f1_size_u > ${WD}/breakseq2/${ID}/${ID}.runfilter.log
				[ ! -d "${WD}/breakseq2" ] && mkdir -p -m 775 "${WD}/breakseq2"
				[ -d "${WD}/breakseq2/${ID}" ] && $PYTHON $BIN/VCFfilterSize.py $breakseq2input $f1_size_l $f1_size_u > ${WD}/breakseq2/${ID}/${ID}.runfilter.log || $PYTHON $BIN/VCFfilterSize.py $breakseq2input $f1_size_l $f1_size_u > ${WD}/breakseq2/${ID}.runfilter.log
				t=$(basename $breakseq2input .gz); breakseq2input="$(dirname $breakseq2input)/${t%.*}${filterName}.vcf"
			fi
			MetaSV_arg="$MetaSV_arg --breakseq_vcf $breakseq2input"
		fi
		if  echo "${MSts[*]}" | grep -w "pindel" &>/dev/null && [ "${premetasv_flag[2]}" = 'TRUE' ]; then
			pindelinput=${pindelvcf:-${WD}/pindel/${ID}/${ID}.pindel.chr${CHR}.vcf}
			if [ "$Raw" = 'FALSE' ]; then
				filterName=".size_${f1_size_l}_${f1_size_u}"
#				$PYTHON $BIN/VCFfilterSize.py $pindelinput $f1_size_l $f1_size_u > ${WD}/pindel/${ID}/${ID}.pindel.chr${CHR}.runfilter.log
				[ ! -d "${WD}/pindel" ] && mkdir -p -m 775 "${WD}/pindel"
				[ -d "${WD}/pindel/${ID}" ] && $PYTHON $BIN/VCFfilterSize.py $pindelinput $f1_size_l $f1_size_u > ${WD}/pindel/${ID}/${ID}.pindel.chr${CHR}.runfilter.log || $PYTHON $BIN/VCFfilterSize.py $pindelinput $f1_size_l $f1_size_u > ${WD}/pindel/${ID}.pindel.chr${CHR}.runfilter.log
				t=$(basename $pindelinput .gz); pindelinput="$(dirname $pindelinput)/${t%.*}${filterName}.vcf"
			fi
			MetaSV_arg="$MetaSV_arg --pindel_vcf $pindelinput"
		fi
		touch ${WD}/metasv/${ID}/${ID}.metasv.chr${CHR}.err
		grep "metasv\.main\s*All Done" ${WD}/metasv/${ID}/${ID}.metasv.chr${CHR}.err &>/dev/null && echo "Seems you have successfully ran metasv before. You may remove ${WD}/${ID}/${ID}.metasv.chr${CHR}.err file, if you are re-submitting for problems regarding metasv results."
		nrr=0
		until grep "metasv\.main\s*All Done" ${WD}/metasv/${ID}/${ID}.metasv.chr${CHR}.err  &>/dev/null || [ "$nrr" -ge "$reRunStop" ]
		do
			if [ -s ${WD}/metasv/${ID}/${ID}.metasv.chr${CHR}.err ]; then
				echo "Error: re-run MetaSV for $ID chr${CHR}!"
				let nrr++
			fi
			$METASV $MetaSV_arg 2>${WD}/metasv/${ID}/${ID}.metasv.chr${CHR}.err
		done
		! grep "metasv\.main\s*All Done" ${WD}/metasv/${ID}/${ID}.metasv.chr${CHR}.err  &>/dev/null && echo "MetaSV failed for chromosome $CHR!" && MetaSVf='FALSE'
	else
		echo "MetaSV: Unknown Chromosome!"; exit 1
	fi

	echo -n "Finish MetaSV @ "; date

	if [ "$Raw" = 'FALSE' ] && [ "$MetaSVf" = 'TRUE' ]; then
		echo -n "Filter for MetaSV results @ "; date
		source /picb/humpopg7/Frq_share/ROOT/bin/thisroot.sh
		echo -n "Performing genotyping for MetaSV output @ "; date
		if [ ! -e "${WD}/cnvnator/${ID}/${ID}.root" ]; then
			echo -n "Seems you don't have ${ID}.root (in ${WD}/cnvnator/${ID}/), I will run CNVnator to make the root file... @ "; date
			[ ! -d ${WD}/cnvnator/${ID}/ ] && mkdir -p -m 775 ${WD}/cnvnator/${ID}/
			$CNVNATOR -root ${WD}/cnvnator/${ID}/${ID}.root -unique -genome GRCh37 -tree $BAM > ${WD}/cnvnator/${ID}/${ID}.cnvnator.runlog
			echo -n "Finish making CNVnator root file @ "; date
		fi
		if [ "$CHR" = 'all' ] || [ "$CHR" = 'autosome' ]; then
			$PYTHON ${BIN}/SeparateFilter_CNVvcf.gamma.py ${WD}/metasv/${ID}/${ID}.MetaSV.vcf ${gap_file} ${gap_frac} ${RD_confSize} -${f2_size_l} +${f2_size_u} > ${WD}/metasv/${ID}/${ID}.MetaSV.runfilter.log
			$PYTHON ${BIN}/MetaSV_selfOverlap_v5.py ${WD}/metasv/${ID}/${ID}.MetaSV.DEL.vcf DEL > ${WD}/metasv/${ID}/${ID}.filtMetaSV.log
			$PYTHON ${BIN}/MetaSV_selfOverlap_v5.py ${WD}/metasv/${ID}/${ID}.MetaSV.DUP.vcf DUP >> ${WD}/metasv/${ID}/${ID}.filtMetaSV.log
			$PERL ${BIN}/Convert_VCF2CNVR.pl ${WD}/metasv/${ID}/${ID}.MetaSV.DEL.nonOverlap.vcf >> ${WD}/metasv/${ID}/${ID}.filtMetaSV.log
			$PERL ${BIN}/Convert_VCF2CNVR.pl ${WD}/metasv/${ID}/${ID}.MetaSV.DUP.nonOverlap.vcf >> ${WD}/metasv/${ID}/${ID}.filtMetaSV.log
			echo "exit" |cat ${WD}/metasv/${ID}/${ID}.MetaSV.DEL.nonOverlap.region.txt - | $CNVNATOR -root ${WD}/cnvnator/${ID}/${ID}.root -genotype $BIN_SIZE 1>${WD}/metasv/${ID}/${ID}.MetaSV.DEL.genotype.txt 2>${WD}/metasv/${ID}/${ID}.filtMetaSV.err 
			echo "exit" |cat ${WD}/metasv/${ID}/${ID}.MetaSV.DUP.nonOverlap.region.txt - | $CNVNATOR -root ${WD}/cnvnator/${ID}/${ID}.root -genotype $BIN_SIZE 1>${WD}/metasv/${ID}/${ID}.MetaSV.DUP.genotype.txt 2>${WD}/metasv/${ID}/${ID}.filtMetaSV.err 
			$PERL ${BIN}/Validate_by_genotype.pl ${WD}/metasv/${ID}/${ID}.MetaSV.DEL.nonOverlap.vcf ${WD}/metasv/${ID}/${ID}.MetaSV.DEL.genotype.txt >> ${WD}/metasv/${ID}/${ID}.filtMetaSV.log 
			$PERL ${BIN}/Validate_by_genotype.pl ${WD}/metasv/${ID}/${ID}.MetaSV.DUP.nonOverlap.vcf ${WD}/metasv/${ID}/${ID}.MetaSV.DUP.genotype.txt >> ${WD}/metasv/${ID}/${ID}.filtMetaSV.log 
			$PYTHON ${BIN}/MetaSV.DEL_DUP.merge.py ${WD}/metasv/${ID}/${ID}.MetaSV.DEL.nonOverlap.reGenoPass.vcf ${WD}/metasv/${ID}/${ID}.MetaSV.DUP.nonOverlap.reGenoPass.vcf >> ${WD}/metasv/${ID}/${ID}.filtMetaSV.log
		elif [[ "$CHR" =~ ^[1-9]$ || "$CHR" =~ ^1[0-9]$ || "$CHR" =~ ^2[0-2]$ || "$CHR" =~ ^[XY]$ ]]; then
			$PYTHON ${BIN}/SeparateFilter_CNVvcf.gamma.py ${WD}/metasv/${ID}/chr${CHR}/variants.vcf.gz ${gap_file} ${gap_frac} ${RD_confSize} -${f2_size_l} +${f2_size_u} > ${WD}/metasv/${ID}/chr${CHR}/MetaSV.runfilter.log
			$PYTHON ${BIN}/MetaSV_selfOverlap_v5.py ${WD}/metasv/${ID}/chr${CHR}/variants.DEL.vcf DEL > ${WD}/metasv/${ID}/chr${CHR}/filtMetaSV.log
			$PYTHON ${BIN}/MetaSV_selfOverlap_v5.py ${WD}/metasv/${ID}/chr${CHR}/variants.DUP.vcf DUP >> ${WD}/meatsv/${ID}/chr${CHR}/filtMetaSV.log
			$PERL ${BIN}/Convert_VCF2CNVR.pl ${WD}/metasv/${ID}/chr${CHR}/variants.DEL.nonOverlap.vcf >> ${WD}/metasv/${ID}/chr${CHR}/filtMetaSV.log
			$PERL ${BIN}/Convert_VCF2CNVR.pl ${WD}/metasv/${ID}/chr${CHR}/variants.DUP.nonOverlap.vcf >> ${WD}/metasv/${ID}/chr${CHR}/filtMetaSV.log
			echo "exit" |cat ${WD}/metasv/${ID}/chr${CHR}/variants.DEL.nonOverlap.region.txt - | $CNVNATOR -root ${WD}/cnvnator/${ID}/${ID}.root -genotype $BIN_SIZE 1>${WD}/metasv/${ID}/chr${CHR}/variants.DEL.genotype.txt 2>${WD}/metasv/${ID}/chr${CHR}/filtMetaSV.err 
			echo "exit" |cat ${WD}/metasv/${ID}/chr${CHR}/variants.DUP.nonOverlap.region.txt - | $CNVNATOR -root ${WD}/cnvnator/${ID}/${ID}.root -genotype $BIN_SIZE 1>${WD}/metasv/${ID}/chr${CHR}/variants.DUP.genotype.txt 2>${WD}/metasv/${ID}/chr${CHR}/filtMetaSV.err 
			$PERL ${BIN}/Validate_by_genotype.pl ${WD}/metasv/${ID}/chr${CHR}/variants.DEL.nonOverlap.vcf ${WD}/metasv/${ID}/chr${CHR}/variants.DEL.genotype.txt >> ${WD}/metasv/${ID}/chr${CHR}/filtMetaSV.log 
			$PERL ${BIN}/Validate_by_genotype.pl ${WD}/metasv/${ID}/chr${CHR}/variants.DUP.nonOverlap.vcf ${WD}/metasv/${ID}/chr${CHR}/variants.DUP.genotype.txt >> ${WD}/metasv/${ID}/chr${CHR}/filtMetaSV.log 
			$PYTHON ${BIN}/MetaSV.DEL_DUP.merge.py ${WD}/metasv/${ID}/chr${CHR}/variants.DEL.nonOverlap.reGenoPass.vcf ${WD}/metasv/${ID}/chr${CHR}/variants.DUP.nonOverlap.reGenoPass.vcf >> ${WD}/metasv/${ID}/chr${CHR}/filtMetaSV.log
		else
			echo "MetaSV: Unknown Chromosome!"; exit 1
		fi
		echo -n "Finish filtering for MetaSV results @ "; date
	fi

	echo -n "Finish entire calling pipeline @ "; date
fi

if [ "$geno" = 'TRUE' ] && [ "$MetaSVf" = 'TRUE' ]; then
	source /picb/humpopg7/Frq_share/ROOT/bin/thisroot.sh
	echo -n "Performing genotyping for MetaSV output @ "; date
	if [ ! -e "${WD}/cnvnator/${ID}/${ID}.root" ]; then
		echo -n "Seems you don't have ${ID}.root (in ${WD}/cnvnator/${ID}/), I will run CNVnator to make the root file... @ "; date
		[ ! -d ${WD}/cnvnator/${ID}/ ] && mkdir -p -m 775 ${WD}/cnvnator/${ID}/
		$CNVNATOR -root ${WD}/cnvnator/${ID}/${ID}.root -unique -genome GRCh37 -tree $BAM > ${WD}/cnvnator/${ID}/${ID}.cnvnator.runlog
		echo -n "Finish making CNVnator root file @ "; date
	fi
		
	if [ "$Raw" = 'FALSE' ]; then
		if [ "$CHR" = 'all' ] || [ "$CHR" = 'autosome' ]; then
			$PERL ${BIN}/Convert_VCF2CNVR.pl ${WD}/metasv/${ID}/${ID}.MetaSV.DEL_DUP.nonOverlap.reGenoPass.vcf >> ${WD}/metasv/${ID}/${ID}.filtMetaSV.log
			echo "exit" |cat ${WD}/metasv/${ID}/${ID}.MetaSV.DEL_DUP.nonOverlap.reGenoPass.region.txt - | $CNVNATOR -root ${WD}/cnvnator/${ID}/${ID}.root -genotype $BIN_SIZE 1>${WD}/metasv/${ID}/${ID}.MetaSV.DEL_DUP.genotype.txt 2>${WD}/metasv/${ID}/${ID}.filtMetaSV.err
		elif [[ "$CHR" =~ ^[1-9]$ || "$CHR" =~ ^1[0-9]$ || "$CHR" =~ ^2[0-2]$ || "$CHR" =~ ^[XY]$ ]]; then
			$PERL ${BIN}/Convert_VCF2CNVR.pl ${WD}/metasv/${ID}/chr${CHR}/variants.DEL_DUP.nonOverlap.reGenoPass.vcf >> ${WD}/metasv/${ID}/chr${CHR}/filtMetaSV.log
			echo "exit" |cat ${WD}/metasv/${ID}/chr${CHR}/variants.DEL_DUP.nonOverlap.reGenoPass.region.txt - | $CNVNATOR -root ${WD}/cnvnator/${ID}/${ID}.root -genotype $BIN_SIZE 1>${WD}/metasv/${ID}/chr${CHR}/variants.DEL_DUP.genotype.txt 2>${WD}/metasv/${ID}/chr${CHR}/filtMetaSV.err
		else
			echo "postMetaSV genotyping: Unknown Chromosome!"; exit 1
		fi
	else
		if [ "$CHR" = 'all' ] || [ "$CHR" = 'autosome' ]; then
			$PERL ${BIN}/Convert_VCF2CNVR.pl ${WD}/metasv/${ID}/${ID}.MetaSV.vcf > ${WD}/metasv/${ID}/${ID}.postMetaSV.log
			echo "exit" |cat ${WD}/metasv/${ID}/${ID}.MetaSV.region.txt - | $CNVNATOR -root ${WD}/cnvnator/${ID}/${ID}.root -genotype $BIN_SIZE 1>${WD}/metasv/${ID}/${ID}.MetaSV.genotype.txt 2>${WD}/metasv/${ID}/${ID}.postMetaSV.err
		elif [[ "$CHR" =~ ^[1-9]$ || "$CHR" =~ ^1[0-9]$ || "$CHR" =~ ^2[0-2]$ || "$CHR" =~ ^[XY]$ ]]; then
			$PERL ${BIN}/Convert_VCF2CNVR.pl ${WD}/metasv/${ID}/chr${CHR}/variants.vcf.gz > ${WD}/metasv/${ID}/chr${CHR}/postMetaSV.log
			echo "exit" |cat ${WD}/metasv/${ID}/chr${CHR}/variants.region.txt - | $CNVNATOR -root ${WD}/cnvnator/${ID}/${ID}.root -genotype $BIN_SIZE 1>${WD}/metasv/${ID}/chr${CHR}/variants.genotype.txt 2>${WD}/metasv/${ID}/chr${CHR}/postMetaSV.err
		else
			echo "postMetaSV genotyping: Unknown Chromosome!"; exit 1
		fi
	fi
	echo -n "Finish genotyping for MetaSV @ "; date
fi

echo -n "Finish entire pipeline for $ID @ "; date

