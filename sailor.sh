#!/bin/sh
module load SAMtools
module load BEDTools
module load Anaconda2/4.0.0
module load BCFtools
for i in `ls *_pass2Aligned.sortedByCoord.out.bam`
do 
bam=~/Aphid_work/45rna/${i}
i=${i/_pass2Aligned.sortedByCoord.out.bam/}

mkdir ${i}_sailorgt
cd ${i}_sailorgt
cp ../example.yaml run.yaml
sed -i "s#bath.bam#${bam}#g"  run.yaml

bsub -J blast -n 1 -R "span[hosts=1] rusage[mem=10GB]" -o %J.out -e %J.err -q q2680v2 "cwltool ~/Aphid_work/sailor-1.2.0/CWL-SINGULARITY-pipeline-building-code/cwl/wf_rnaediting2strands.cwl  run.yaml"
cd ..

done
