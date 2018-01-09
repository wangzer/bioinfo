# SNP Calling Using samtools
#! /bin/bash
set -u
set -e
set -o pipefail
# set work path
# PATH=/bin:/usr/bin:/usr/sbin:/usr/local/bin:/usr/local/sbin:~/bin
# or if you have install miniconda/anaconda and configure the bioconda 
source activa biostar
work_dir=$1
sample_info=$2
reference=${work_dir}/data/reference/TAIR10_chr_all.fas
# alignment
bwa index $reference
sample_names=($(cut -f 1  "$sample_info" | uniq))
mkdir -p ${work_dir}/analysis/align
cd ${work_dir}/data/seq
for sample in ${sample_names[@]}
do
    # create an output file from the sample name
    result_file="${sample}.sam"
    bwa mem $reference ${sample}_R1.fastq ${sample}_R2.fastq > ${work_dir}/analysis/align/$result_file
    #echo "$result_file"
done
#snp calling
cd ${work_dir}/analysis/align
samfiles=$(ls *sam)
echo $samfiles
mkdir -p ${work_dir}/analysis/snp
for file in ${samfiles[@]}
do
    output=${file%.*}
    samtools view -b -o ${output}.bam ${file}
    samtools sort -o ${output}.sorted.bam ${output}.bam
    samtools index ${output}.sorted.bam
    samtools mpileup -uv -t AD,DP -f $reference ${output}.sorted.bam | bcftools call -vm -Ov > ../snp/${output}.vcf
done
# convert vcf4.2 to vcf 4.1
cd ${work_dir}/analysis/snp
infiles=$(ls *vcf)
for infile in ${infiles[@]}
do
    bcftools  view --max-alleles 2 -O v ${infile} | \
    sed "s/##fileformat=VCFv4.2/##fileformat=VCFv4.1/" | \
    sed "s/(//" | \
    sed "s/)//" | \
    sed "s/,Version=\"3\">/>/" | \
    bcftools view -O z > ${infile%%.*}.dg.vcf.gz
done
