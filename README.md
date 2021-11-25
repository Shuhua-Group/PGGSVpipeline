# PGGSVpipeline

PGGSVpipeline is an automatic workflow for NGS SV calling.

## PGGSVpipeline usage

```
usage:   PGGSVpipeline.sh [options]

PGG.SVpipeline.sh --bam(f)ile indi.realigned.bam [[ --sample(i)d sampleid --(c)hr chromosome --work(d)ir working dir --(s)teps steps to run the pipeline --(m)etasvtools tools used to merge SVs by MetaSV --(r)eadlength readlength --(a)ssembly trun on assembly for MetaSV --(R)aw output raw results, instead of an embedded fitlering ] --dry specify running tools in dry-run model(, for tools already DONE) --(F)orce force to run MetaSV if any pre-metasv tool fails --(G)enotype perform genotyping after MetaSV step using CNVnator --gap_file genomic gaps to filter the CNV regions --gap_frac cutoff of gap filtering --pindel_vcf pindel vcf for MetaSV --breakdancer_naive breakdancer results for MetaSV --cnvnator_naive cnvnator results for MetaSV --breakseq2_vcf breakseq2 vcf for MetaSV --lumpy_vcf lumpy vcf for MetaSV ]
```
