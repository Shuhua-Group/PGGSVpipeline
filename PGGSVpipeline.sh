#!/bin/bash
SD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" #Script_Directory
source ${SD}/parameter.sh
WD=`pwd` #Work_Directory

function usage0(){
	cat << INTRO

PGGSVpipeline version 2.0, 20220826.

Program:   PGGSVpipeline (Population Genomics Group Structural Variation - calling pipeline)
Contact:   Yimin Wang <wangyimin@picb.ac.cn>

Usage:     PGGSVpipeline.sh -b <input.bam> -i <sampleSID> 
           -s <software_list> -r <reference_genome> -f <ref.fasta>

Required parameters:
           -b    the mapped sequences in bam format
           -r    choose a build-in reference genome (hg19/hg38/chm13)
        or
           -f    use a user reference genome


Optional parameters:
           -i    task SID
           -s    enter a comma-separated list of software. 'all' is equal to 'breakdancer,breakseq2,cnvnator,lumpy,manta,metasv'
                 please note that breakseq2 can only be used with hg19 or hg38 (default all)

INTRO
}

######################## functions #########################
function run_breakdancer(){
	BAM=$1; ID=$2; WD=$3
	${BREAKDANCER}/perl ${BREAKDANCER}/bam2cfg.pl ${BAM} > ${WD}/${ID}.breakdancer.cfg

	for chr in ${CHRlist[@]}
	do
		${BREAKDANCER}/breakdancer-max -o ${chr} ${WD}/${ID}.breakdancer.cfg > ${WD}/${ID}.${chr}.SV.output &
	done
	wait
	CHRbox=($CHRlist)
	if [ -s ${WD}/${ID}.breakdancer.cfg -a -s ${WD}/${ID}.${CHRbox[20]}.SV.output ]; then
		echo -n "${ID} breakdancer DONE @ "; date
		echo "${ID} breakdancer DONE @ " > ${WD}/taskdone.log
		date >> ${WD}/taskdone.log
	else
		echo -n "${ID} breakdancer FAILED @ "; date
		touch ${WD}/${ID}.breakdancer.fail
	fi
}

function run_breakseq2(){
	BAM=$1; ID=$2; WD=$3
	touch ${WD}/${ID}.breakseq2.runlog
	nrr=0
	until grep "breakseq2_workflow   Done" ${WD}/${ID}.breakseq2.runlog &>/dev/null || [ "${nrr}" -ge "${reRunStop}" ]
	do
		if [ ${nrr} -gt 0 ]; then
			echo "Error: re-run breakseq2 for ${ID}!"
		fi
		${BREAKSEQ2}/python ${BREAKSEQ2}/run_breakseq2.py --bams ${BAM} --work ${WD} --reference ${USER_REF} --chromosomes ${CHRlist} --sample ${ID} --bplib_gff ${BS2_bplib} --bwa ${BREAKSEQ2}/bwa --samtools ${BREAKSEQ2}/samtools --min_span 10 --window ${BIN_SIZE} --min_overlap 10 1>/dev/null 2>${WD}/${ID}.breakseq2.runlog
		let nrr++
	done
	if grep "breakseq2_workflow   Done" ${WD}/${ID}.breakseq2.runlog &>/dev/null; then
		echo -n "${ID} breakseq2 DONE @ "; date
		echo "${ID} breakseq2 DONE @ " > ${WD}/taskdone.log
		date >> ${WD}/taskdone.log
	else
		echo -n "${ID} breakseq2 FAILED @ "; date
		touch ${WD}/${ID}.breakseq2.fail
	fi
}

function run_cnvnator(){
	export PATH=${CNVNATOR}:$PATH
	BAM=$1; ID=$2; WD=$3
	touch ${WD}/${ID}.cnvnator.runlog
	if [ "${DEF_REF}" = 'hg19' ]; then
		${CNVNATOR}/cnvnator -root ${WD}/${ID}.root -genome hg19 -chrom ${CHRlist} -tree ${BAM} > ${WD}/${ID}.cnvnator.runlog
		${CNVNATOR}/cnvnator -root ${WD}/${ID}.root -genome hg19 -chrom ${CHRlist} -his ${BIN_SIZE} -fasta ${USER_REF} > ${WD}/${ID}.cnvnator.runlog
	else
		${CNVNATOR}/cnvnator -root ${WD}/${ID}.root -tree ${BAM} > ${WD}/${ID}.cnvnator.runlog
		${CNVNATOR}/cnvnator -root ${WD}/${ID}.root -his ${BIN_SIZE} -chrom ${CHRlist} -fasta ${USER_REF} > ${WD}/${ID}.cnvnator.runlog
	fi
	${CNVNATOR}/cnvnator -root ${WD}/${ID}.root -stat ${BIN_SIZE} > ${WD}/${ID}.cnvnator.runlog
	${CNVNATOR}/cnvnator -root ${WD}/${ID}.root -partition ${BIN_SIZE} > ${WD}/${ID}.cnvnator.runlog
	${CNVNATOR}/cnvnator -root ${WD}/${ID}.root -call ${BIN_SIZE} > ${WD}/${ID}.rawcnv
	if [ -s ${WD}/${ID}.rawcnv ]; then
		echo -n "${ID} cnvnator DONE @ "; date
		echo "${ID} cnvnator DONE @ " > ${WD}/taskdone.log
		date >> ${WD}/taskdone.log
	else
		echo -n "${ID} cnvnator FAILED @ "; date
		touch ${WD}/${ID}.cnvnator.fail
	fi
}

function run_lumpy(){
	export PATH=${LUMPY}:$PATH
	BAM=$1; ID=$2; WD=$3
	${SAMTOOLS} view -b -F 1294 ${BAM} > ${WD}/${ID}.discordants.unsorted.bam &
	${SAMTOOLS} view -h ${BAM} | ${LUMPY}/extractSplitReads_BwaMem -i stdin | ${SAMTOOLS} view -Sb - > ${WD}/${ID}.splitters.unsorted.bam &
	wait
	${SAMTOOLS} sort ${WD}/${ID}.discordants.unsorted.bam -o ${WD}/${ID}.discordants.bam &
	${SAMTOOLS} sort ${WD}/${ID}.splitters.unsorted.bam -o ${WD}/${ID}.splitters.bam &
	wait
	touch ${WD}/${ID}.lumpy.runlog
	nrr=0
	until grep "LUMPY Express done" ${WD}/${ID}.lumpy.runlog &>/dev/null || [ "${nrr}" -ge "${reRunStop}" ]
	do
		if [ ${nrr} -gt 0 ]; then
			echo "Error: re-run Lumpy for ${ID}!"
		fi
		${LUMPY}/lumpyexpress -B ${BAM} -S ${WD}/${ID}.splitters.bam -D ${WD}/${ID}.discordants.bam -o ${WD}/${ID}.lumpy.vcf > ${WD}/${ID}.lumpy.runlog
		let nrr++
	done
	rm ${WD}/${ID}.discordants.unsorted.bam ${WD}/${ID}.splitters.unsorted.bam
	if grep "LUMPY Express done" ${WD}/${ID}.lumpy.runlog &>/dev/null; then
		echo -n "${ID} lumpy DONE @ "; date
		echo "${ID} lumpy DONE @ " > ${WD}/taskdone.log
		date >> ${WD}/taskdone.log
	else
		echo  -n "${ID} lumpy FAILED @ "; date
		touch ${WD}/${ID}.lumpy.fail
	fi
}

function run_manta(){
	BAM=$1; ID=$2; WD=$3
	${MANTA}/python ${MANTA}/configManta.py --bam ${BAM} --referenceFasta ${USER_REF} --runDir ${WD}/ &>/dev/null
	${MANTA}/python ${WD}/runWorkflow.py -j 10 &>/dev/null
	if [ -s ${WD}/results/variants/diploidSV.vcf.gz ]; then
		echo -n "${ID} manta DONE @ "; date
		echo "${ID} manta DONE @ " > ${WD}/taskdone.log
		date >> ${WD}/taskdone.log
	else
		echo -n "${ID} manta FAILED @ "; date
		touch ${WD}/${ID}.manta.fail
	fi
}

######################## parameters #########################
if [ x$1 != x ]
then
	while getopts "b:i:s:r:f:" arg
	do
		case $arg in
		b)
			BAM=${OPTARG}
			;;
		i)
			SID=${OPTARG}
			;;
		s)
			SOFT=${OPTARG}
			;;
		r)
			DEF_REF=${OPTARG}
			;;
		f)
			USER_REF=${OPTARG}
			;;
		?)
			echo "FAIL: Unknown argument. Please check again."
			exit 1
			;;
		esac
	done
	# check necessary parameters
	if [ -z "${BAM}" ] || [ ! -e "${BAM}" ];then
		echo "FAIL: Not available BAM file ${BAM}!"
		echo "####################################"
		usage0; exit 1
	fi

	if [ -z "${SID}" ]; then
		SID=$(basename ${BAM} .bam)
		SID=$(basename ${SID} .realigned)
		SID=$(basename ${SID} .dedup)
		SID=$(basename ${SID} .bqsr)
	fi

	if [ -z "${DEF_REF}" ] && [ -z "${USER_REF}" ]; then
		echo "FAIL: REF not selected!"
		echo "####################################"
		usage0; exit 1
	fi

	if [ ! -z "${USER_REF}" ] && [ ! -e "${USER_REF}" ]; then
		echo "FAIL: Not available REF ${USER_REF}!"
		echo "####################################"
		usage0; exit 1
	fi

	if [ ! -z "${USER_REF}" ]; then
		[ ${SOFT} = "all" ] && SOFT="breakdancer,cnvnator,lumpy,manta,metasv"
	else
		if [ "${DEF_REF}" = 'hg19' ]; then
			USER_REF="${SD}/resources/human_g1k_v37.fasta"
			BS2_bplib="${SD}/resources/bplib.hg19.gff"
			CHRlist=${HG19_CHRlist}
			if [ "${SOFT}" = 'all' ]; then
				SOFT="breakdancer,breakseq2,cnvnator,lumpy,manta,metasv"
			fi
		elif [ "${DEF_REF}" = 'hg38' ]; then
			USER_REF="${SD}/resources/GRCh38_full_analysis_set_plus_decoy_hla.fa"
			BS2_bplib="${SD}/resources/bplib.hg38.gff"
			if [ "${SOFT}" = 'all' ]; then
				SOFT="breakdancer,breakseq2,cnvnator,lumpy,manta,metasv"
			fi
		elif [ "${DEF_REF}" = 'chm13' ]; then
			USER_REF="${SD}/resources/chm13v2.0.fa"
			if [ "${SOFT}" = 'all' ]; then
				SOFT="breakdancer,cnvnator,lumpy,manta,metasv"
			fi
		else
			echo "FAIL: Not available REF ${DEF_REF}!"
			echo "####################################"
			usage0; exit 1
		fi
	fi

	#rebuild bam index
	mkdir -p -m 775 "${WD}/bam/" 2>/dev/null
	if [ ! -e "${WD}/bam/${SID}.bam.bai" ]; then
		ln -s ${BAM} ${WD}/bam/${SID}.bam
		cd ${WD}/bam/
		${SAMTOOLS} index -b ${SID}.bam
		cd ${WD}
	fi
	BAM="${WD}/bam/${SID}.bam"

	# get softs
	SOFTS=(${SOFT//,/ })
	run=()
	for step in ${SOFTS[@]}
	do
		if echo "${tools[@]}" |grep -w "${step}" &>/dev/null; then
			run+=(${step})
		fi
	done
	run_MetaSV='FALSE'
	if echo "${SOFTS[@]}" |grep -w "metasv" &>/dev/null; then
		run_Metasv='TRUE'
	fi

	#running
	for tool in ${run[@]}
	do
		if [ "$tool" = 'breakdancer' ]; then
			if [ -e "${WD}/breakdancer/${SID}/taskdone.log" ]; then
				echo -n "${SID} breakdancer EXIST. Skip running... @ "; date
			else
				[ -d "${WD}/breakdancer/${SID}" ] && rm -r ${WD}/breakdancer/${SID}
				mkdir -p -m 775 "${WD}/breakdancer/${SID}"
				run_breakdancer ${BAM} ${SID} ${WD}/breakdancer/${SID} 2>${WD}/${SID}.breakdancer.err &
				echo -n "${SID} running breakdancer... @ "; date
			fi
		elif [ "$tool" = 'breakseq2' ]; then
			if [ -e "${WD}/breakseq2/${SID}/taskdone.log" ]; then
				echo -n "${SID} breakseq2   EXIST. Skip running... @ "; date
			else
				[ -d "${WD}/breakseq2/${SID}" ] && rm -r ${WD}/breakseq2/${SID}
				mkdir -p -m 775 "${WD}/breakseq2/${SID}"
				run_breakseq2 ${BAM} ${SID} ${WD}/breakseq2/${SID} 2>${WD}/${SID}.breakseq2.err &
				echo -n "${SID} running breakseq2... @ "; date
			fi
		elif [ "$tool" = 'cnvnator' ]; then
			if [ -e "${WD}/cnvnator/${SID}/taskdone.log" ]; then
				echo -n "${SID} cnvnator    EXIST. Skip running... @ "; date
			else
				[ -d "${WD}/cnvnator/${SID}" ] && rm -r ${WD}/cnvnator/${SID}
				mkdir -p -m 775 "${WD}/cnvnator/${SID}" 
				run_cnvnator ${BAM} ${SID} ${WD}/cnvnator/${SID} 2>${WD}/${SID}.cnvnator.err &
				echo -n "${SID} running cnvnator... @ "; date
			fi
		elif [ "$tool" = 'lumpy' ]; then
			if [ -e "${WD}/lumpy/${SID}/taskdone.log" ]; then
				echo -n "${SID} lumpy       EXIST. Skip running... @ "; date
			else
				[ -d "${WD}/lumpy/${SID}" ] && rm -r ${WD}/lumpy/${SID}
				mkdir -p -m 775 "${WD}/lumpy/${SID}" 
				run_lumpy ${BAM} ${SID} ${WD}/lumpy/${SID} 2>${WD}/${SID}.lumpy.err &
				echo -n "${SID} running lumpy... @ "; date
			fi
		elif [ "$tool" = 'manta' ]; then
			if [ -e "${WD}/manta/${SID}/taskdone.log" ]; then
				echo -n "${SID} manta       EXIST. Skip running... @ "; date
			else
				[ -d "${WD}/manta/${SID}" ] && rm -r ${WD}/manta/${SID}
				mkdir -p -m 775 "${WD}/manta/${SID}"
				run_manta ${BAM} ${SID} ${WD}/manta/${SID} 2>${WD}/${SID}.manta.err &
				echo -n "${SID} running manta... @ "; date
			fi
		fi
	done
	wait

	## Check failtures
	result_check=(TRUE TRUE TRUE TRUE TRUE)
	for tool in ${MSts[@]}
	do
		[ "$tool" = "breakdancer" ] && [ -e "${WD}/breakdancer/${SID}/${SID}.breakdancer.fail" ] && result_check[0]='FALSE'
		[ "$tool" = "breakseq2" ] && [ -e "${WD}/breakseq2/${SID}/${SID}.breakseq2.fail" ] && result_check[1]='FALSE'
		[ "$tool" = "cnvnator" ] && [ -e "${WD}/cnvnator/${SID}/${SID}.cnvnator.fail" ] && result_check[2]='FALSE'
		[ "$tool" = "lumpy" ] && [ -e "${WD}/lumpy/${SID}/${SID}.lumpy.fail" ] && result_check[3]='FALSE'
		[ "$tool" = "manta" ] && [ -e "${WD}/manta/${SID}/${SID}.manta.fail" ] && result_check[4]='FALSE'
	done

	## run embedded filter process and merge by MetaSV
	if echo "${result_check[@]}" |grep -w "FALSE" &>/dev/null; then
		echo "Error: Stop MetaSV merging, due to some tools failing to run successfully, check above logs for help!"
		echo -n "Abandon pipeline @ "; date
	else
		echo -n "Start running MetaSV @ "; date
		CHRs=$(echo {1..22} X Y)
		[ -d "${WD}/metasv/${SID}" ] && rm -r ${WD}/metasv/${SID}
		mkdir -p -m 775 "${WD}/metasv/${SID}"
		MetaSV_arg_ini="--reference ${USER_REF} --sample ${SID} --disable_assembly --num_threads 1 --enable_per_tool_output --keep_standard_contigs --mean_read_length ${readlength}"
		for chr in ${CHRlist[@]}
		do
			MetaSV_arg="${MetaSV_arg_ini} --outdir ${WD}/metasv/${SID}/${chr} --workdir ${WD}/metasv/${SID}/${chr}/tmp_work --chromosomes ${chr}"
			#breakdancer
			if echo "${run[@]}" | grep -w "breakdancer" &>/dev/null; then
				breakdancerinput="${WD}/breakdancer/${SID}/${SID}.${chr}.SV.output"
				${SD}/PGGSV.010.SizeFilter_BreakDancer.sh ${breakdancerinput} ${f1_size_l} ${f1_size_u} > ${WD}/breakdancer/${SID}/${SID}.${chr}.runfilter.log
				breakdancerinput="${WD}/breakdancer/${SID}/${SID}.${chr}.size_${f1_size_l}_${f1_size_u}.SV.output"
				MetaSV_arg="${MetaSV_arg} --breakdancer_native ${breakdancerinput}"
			fi
			#breakseq2
			if echo "${run[@]}" | grep -w "breakseq2" &>/dev/null; then
				breakseq2input=${WD}/breakseq2/${SID}/breakseq.vcf.gz
				${PYTHON} ${SD}/PGGSV.011.SizeFilter.py ${breakseq2input} ${f1_size_l} ${f1_size_u} > ${WD}/breakseq2/${SID}/${SID}.runfilter.log
				breakseq2input="${WD}/breakseq2/${SID}/${SID}.breakseq.size_${f1_size_l}_${f1_size_u}.vcf"
				MetaSV_arg="${MetaSV_arg} --breakseq_vcf ${breakseq2input}"
			fi
			#cnvnator
			if echo "${run[@]}" | grep -w "cnvnator" &>/dev/null; then
				cnvnatorinput="${WD}/cnvnator/${SID}/${SID}.rawcnv"
				MetaSV_arg="${MetaSV_arg} --cnvnator_native ${cnvnatorinput}"
			fi
			#lumpy
			if echo "${run[@]}" | grep -w "lumpy" &>/dev/null; then
				lumpyinput=${WD}/lumpy/${SID}/${SID}.lumpy.vcf
				filterName=".size_${f1_size_l}_${f1_size_u}"
				${PYTHON} ${SD}/PGGSV.011.SizeFilter.py ${lumpyinput} ${f1_size_l} ${f1_size_u} > ${WD}/lumpy/${SID}/${SID}.lumpy.runfilter.log
				lumpyinput="${WD}/lumpy/${SID}/${SID}.lumpy.size_${f1_size_l}_${f1_size_u}.vcf"
				MetaSV_arg="${MetaSV_arg} --lumpy_vcf ${lumpyinput}"
			fi

			#run metasv
			(nrr=0
			until grep "metasv\.main\s*All Done" ${WD}/metasv/${SID}/${SID}.metasv.${chr}.log &>/dev/null || [ "${nrr}" -ge "${reRunStop}" ]
			do
				if [ ${nrr} -gt 0 ]; then
					echo "Error: re-run MetaSV for ${SID} ${chr}!"
					let nrr++
				fi
				${METASV}/python ${METASV}/run_metasv.py ${MetaSV_arg} &>${WD}/metasv/${SID}/${SID}.metasv.${chr}.log
			done) &
		done
		wait

		FAIL='FALSE'
		for chr in ${CHRlist[@]}
		do
			grep "metasv\.main\s*All Done" ${WD}/metasv/${SID}/${SID}.metasv.${chr}.log &>/dev/null || FAIL='TRUE'
		done

		echo -n "Finish MetaSV @ "; date
		if [ "${FAIL}" = 'FALSE' ]; then
			${SD}/PGGSV.020.MergeChr.pl ${WD}/metasv/${SID} ${SID} > ${WD}/metasv/${SID}/${SID}.merge.log
		else
			echo "MetaSV failed for some chromosomes."
		fi
	fi

	#merge
	[ -d "${WD}/CombineOutput/${SID}" ] && rm -r ${WD}/CombineOutput/${SID}
	mkdir -p -m 775 "${WD}/CombineOutput/${SID}"

	${PYTHON} ${SD}/PGGSV.021.MergeManta.py ${WD}/manta/${SID}/results/variants/diploidSV.vcf.gz ${WD}/metasv/${SID}/${SID}.MetaSV.vcf ${WD}/CombineOutput/${SID}/${SID}.MetaSV.vcf
	${PYTHON} ${SD}/PGGSV.022.SeparateFilter.py ${WD}/CombineOutput/${SID}/${SID}.MetaSV.vcf NULL 0 ${RD_confSize} -${f2_size_l} +${f2_size_u} > ${WD}/CombineOutput/${SID}/${SID}.MetaSV.runfilter.log
	${PYTHON} ${SD}/PGGSV.023.SelfOverlap.py ${WD}/CombineOutput/${SID}/${SID}.MetaSV.DEL.vcf DEL > ${WD}/CombineOutput/${SID}/${SID}.filtMetaSV.log
	${PYTHON} ${SD}/PGGSV.023.SelfOverlap.py ${WD}/CombineOutput/${SID}/${SID}.MetaSV.DUP.vcf DUP >> ${WD}/CombineOutput/${SID}/${SID}.filtMetaSV.log
	${PERL} ${SD}/PGGSV.024.vcf2cnvr.pl ${WD}/CombineOutput/${SID}/${SID}.MetaSV.DEL.nonOverlap.vcf >> ${WD}/CombineOutput/${SID}/${SID}.filtMetaSV.log
	${PERL} ${SD}/PGGSV.024.vcf2cnvr.pl ${WD}/CombineOutput/${SID}/${SID}.MetaSV.DUP.nonOverlap.vcf >> ${WD}/CombineOutput/${SID}/${SID}.filtMetaSV.log
	${PYTHON} ${SD}/PGGSV.025.MergeCNV.py ${WD}/CombineOutput/${SID}/${SID}.MetaSV.DEL.nonOverlap.vcf ${WD}/CombineOutput/${SID}/${SID}.MetaSV.DUP.nonOverlap.vcf >> ${WD}/CombineOutput/${SID}/${SID}.filtMetaSV.log

	if [ -e "${WD}/CombineOutput/${SID}/${SID}.MetaSV.INV.vcf" ]; then
		${PYTHON} ${SD}/PGGSV.026.CombineOutput.py ${WD}/CombineOutput/${SID}/${SID}.MetaSV.DEL_DUP.nonOverlap.vcf ${WD}/CombineOutput/${SID}/${SID}.MetaSV.INV.vcf ${WD}/manta/${SID}/results/variants/diploidSV.vcf.gz ${WD}/CombineOutput/${SID}/${SID}.MetaSV.SV.vcf
	else
		${PYTHON} ${SD}/PGGSV.026.CombineOutput_noINV.py ${WD}/CombineOutput/${SID}/${SID}.MetaSV.DEL_DUP.nonOverlap.vcf ${WD}/manta/${SID}/results/variants/diploidSV.vcf.gz ${WD}/CombineOutput/${SID}/${SID}.MetaSV.SV.vcf
	fi
	
	#clean
	PGGSV.030.DocCleaner.sh -s ${SID}
	echo -n "Finish entire calling pipeline @ "; date
else
	usage0
fi