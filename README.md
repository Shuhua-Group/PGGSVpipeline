# PGGSVpipeline

PGGSVpipeline version 2.0, 20220826.

* Program:   PGGSVpipeline (Population Genomics Group Structural Variation calling pipeline)
* Contact:   Yimin Wang <wangyimin@picb.ac.cn>
```
Usage:     PGGSVpipeline.sh -b <input.bam> -i <sampleSID> 
           -s <software_list> -r <reference_genome> -f <ref.fasta>
```
Required parameters:
```
           -b   the mapped sequences in bam format
           -r   choose a build-in reference genome (hg19/hg38/chm13)
        or
           -f   use a user reference genome
```

Optional parameters:
```
           -i   task SID
           -s   enter a comma-separated list of software. 'all' is equal to 'breakdancer,breakseq2,cnvnator,lumpy,manta,metasv'
                please note that breakseq2 can only be used with hg19 or hg38 (default all)
```

Please replace the path in the parameter.sh file in the script directory with your own path. We recommend using conda to install:
> breakdancer

> breakseq2

> cnvnator

> lumpy

> manta

> metasv

to avoid errors due to conflicting dependent environments.
