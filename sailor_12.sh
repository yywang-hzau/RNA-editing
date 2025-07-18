dna=$1
rna=$2
sample=$3
dnadir='/public/home/yywang/5hosts_2021/'
bed=${dnadir}${dna}.bed
cd  ${sample}
infile=${rna}_deduplicate_split.fwd.sorted.rmdup.readfiltered.formatted.vcf
vcftools --vcf ${infile} --minDP 5 --min-meanDP 5 --recode  --recode-INFO-all --out ${rna}_fwd_filtered
vcf=${rna}_fwd_filtered.recode.vcf
python ../filter_known_snp.py --input ${vcf} --known ${bed} --output ${rna}_fwd_remove_known.vcf 
python ../rank_edits.py -i ${rna}_fwd_remove_known.vcf   -o ${rna}_fwd.conf

infile2=${rna}_deduplicate_split.rev.sorted.rmdup.readfiltered.formatted.vcf
vcftools --vcf ${infile2} --minDP 5 --min-meanDP 5 --recode  --recode-INFO-all --out ${rna}_rev_filtered
vcf2=${rna}_rev_filtered.recode.vcf
python ../filter_known_snp.py --input ${vcf2} --known ${bed} --output ${rna}_rev_remove_known.vcf 
python ../rank_edits.py -i ${rna}_rev_remove_known.vcf   -o ${rna}_rev.conf


