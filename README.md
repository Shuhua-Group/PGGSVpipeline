# PGGSVpipeline

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

## Usage
```
    PGG.SVpipeline.sh -f|--bamfile individual.realigned.bam [options]
```
## options
```
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
    -h 
        print short help
    --help
        print long help
```
## [Note]
```
0.) this pipeline is especially useful to run from BAM files, throughout to get VCF files (in whole genome model). If you have existed pre-metasv results, we strongly suggest you to directly run MetaSV pipeline per se, otherwise we only provide single chromosome model for this situation! Or in a special situation, you may try '--dry' option.
    1.) although we provde 'chr' option, we only support 'whole genome' running for cnvnator, lumpy, and breakseq2.
    2.) breakdancer and pindel are run in divided chromosomes, even if you run in 'whole genome' model, considering the running time
    3.) all the runing parameters for involved tools are as default set, if you are an advanced user, you may modify the running command for any tool in the script, but you should be very clear about what you are doing.
    4.) there are some auto-check for failture running of some tools and will automatically re-run (till to success or server down) -- which is a double-edged sword
    5.) there is a also a stop marker (reRunStop) for 're-run' process, the default time is 5, which means if it continues fail, after 5 times of 're-run'-try, I will not try more, and you may manually check the error and then re-run the whole pipeline and you may consider using '--dry' option.
```
## [Example]
```
    PGG.SVpipeline.sh -f indi.realigned.bam 1>indi.SVpipeline.log 2>indi.SVpipeline.err
```
## [Credit]
    First written by Fu Ruiqing, with contributions from Wu Zhendong, Lou Haiyi, and Zhang Xi. By courtesy of Tian Lei.
    Edited and improved by Wang Yimin & Xie bo
## [Contact]
```
    wangyimin@picb.ac.cn
```
