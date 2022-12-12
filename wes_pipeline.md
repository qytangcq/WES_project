```sh
mkdir wes_project
wd=path_to_work_directory/wes_project
cd $wd
mkdir {raw,clean,qc,align,mutation,bamqc} ##with all fastq file saved in raw fold
```

### Quality Control
- summary stats of raw fastq data
```sh
conda activate wes
cd qc && mkdir {raw,clean}
fastqc $wd/raw/* -t 10 -o $wd/qc/raw   ##10 threads
multiqc -n raw $wd/qc/raw              ##summary reports
```

- trim low quality reads, adapters, etc.
```sh
#!/bin/bash
wd=path_to_work_directory/wes_project
source ~/anaconda3/bin/activate wes
for ID in {1..42}
do

trim_galore --paired -q 20 --phred33 --length 36 -e 0.1 --gzip -j 4 -o $wd/clean $wd/raw/${ID}_*R1* $wd/raw/${ID}_*R2*

done
```

```sh
nohup bash qc.sh &
```

- summary stats of clean fastq data
```sh
conda activate wes
fastqc $wd/clean/* -t 10 -o $wd/qc/clean   ##10 threads
multiqc -n clean $wd/qc/clean              ##summary reports
```

### Align to Reference of Human Genome (hg38)

- buid index for human genome reference
```sh
##reference data required by GATK should be prepared
gatkresource=path_to_gatk_bundle
cd ${gatkresource}
conda activate wes
gzip -d Homo_sapiens_assembly38.fasta.gz
bwa index Homo_sapiens_assembly38.fasta ##creating index file .amb .ann .bwt .pac .sa
```
- Do alignment
```sh
#!/bin/bash
gatkresource=path_to_gatk_bundle
INDEX=${gatkresource}/Homo_sapiens_assembly38.fasta
wd=path_to_work_directory/wes_project
source ~/anaconda3/bin/activate wes
cd $wd/align
ls ../clean/*1.fq.gz | sort > 1
ls ../clean/*2.fq.gz | sort > 2
cut -d"/" -f 3 1 | cut -d"_" -f 1 > 0
paste 0 1 2 > config
cat config | while read id
do
arr=($id)
fq1=${arr[1]}
fq2=${arr[2]}
sample=${arr[0]}
bwa mem -t 5 -R "@RG\tID:${sample}\tSM:${sample}\tLB:WGS\tPL:ILLUMINA" $INDEX $fq1 $fq2 | samtools sort -@ 5 -o ${sample}.bam -
done
```
```
nohup bash align.sh &
```
- summary stats of .bam files
```sh
#!/bin/bash
wd=path_to_work_directory
gatkresource=path_to_gatk_bundle
INDEX=${gatkresource}/Homo_sapiens_assembly38.fasta
cd $wd/align
source ~/anaconda3/bin/activate wes
cat 0 | while read id
do
bam=./${id}.bam
samtools stats -@ 16 --reference ${INDEX} ${bam} > $wd/bamqc/${id}.stat
plot-bamstats -p $wd/bamqc/${id}  $wd/bamqc/${id}.stat
done
```
```
nohup bash bamqc.sh &
```

### GATK Best Practice to Call Germline Variants
```sh
#!/bin/bash
gatkresource=path_to_gatk_bundle
ref=${gatkresource}/Homo_sapiens_assembly38.fasta
snp=${gatkresource}/dbsnp_146.hg38.vcf.gz
indel=${gatkresource}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
wd=path_to_work_directory/wes_project
source ~/anaconda3/bin/activate wes

cd $wd/mutation

find $wd/align -name "*.bam" | while read id; do basename ${id} .bam; done > config

cat config | while read id
do
sample=${id}
gatk --java-options "-Xmx40G -Djava.io.tmpdir=./" MarkDuplicates \
	-I ${sample}.bam  \
	-O ${sample}_marked.bam \
	-M ${sample}.metrics \
	1>${sample}.log.mark 2>&1

gatk --java-options "-Xmx40G -Djava.io.tmpdir=./" FixMateInformation \
	-I ${sample}_marked.bam  \
	-O ${sample}_marked_fixed.bam \
	-SO coordinate  \
	1>${sample}.log.fix 2>&1

samtools index ${sample}_marked_fixed.bam

gatk --java-options "-Xmx40G -Djava.io.tmpdir=./" BaseRecalibrator \
	-R ${ref} \
	-I ${sample}_marked_fixed.bam \
	--known-sites $snp \
	--known-sites $indel \
	-O ${sample}_recal.table \
	1>${sample}.log.recal 2>&1

gatk --java-options "-Xmx40G -Djava.io.tmpdir=./" ApplyBQSR \
	-R ${ref} \
	-I ${sample}_marked_fixed.bam \
	-bqsr ${sample}_recal.table \
	-O ${sample}_bqsr.bam \
	1>${sample}.log.applybqsr 2>&1

gatk --java-options "-Xmx40G -Djava.io.tmpdir=./" HaplotypeCaller \
	-R ${ref} \
	-I ${sample}_bqsr.bam \
	--dbsnp ${snp} \
	-ERC GVCF
	-O ${sample}_raw.g.vcf.gz \
	1>${sample}.log.haplotypecaller 2>&1
done

find ./ -name "*.g.vcf.gz" > gvcf.list

gatk --java-options "-Xmx40G -Djava.io.tmpdir=./" CombineGVCFs \
	-R ${ref} \
	--variant gvcf.list \
	-O combined.g.vcf.gz \
	1>log.combinevcf 2>&1

gatk --java-options "-Xmx40G -Djava.io.tmpdir=./" GenotypeGVCFs \
	-R ${ref} \
	-V combined.g.vcf.gz \
	-O combined.raw.vcf.gz \
	1>log.genotypegvcf 2>&1

## variants filtration with VQSR
gatk --java-options "-Xmx40G -Djava.io.tmpdir=./" VariantRecalibrator \
	-R ${ref} \
	-V combined.raw.vcf.gz \
	-resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${gatkresource}/hapmap_3.3.hg38.vcf.gz \
	-resource:omini,known=false,training=true,truth=false,prior=12.0 ${gatkresource}/1000G_omni2.5.hg38.vcf.gz \
	-resource:1000G,known=false,training=true,truth=false,prior=10.0 ${gatkresource}/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
	-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${gatkresource}/dbsnp_146.hg38.vcf.gz \
	-an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
	-mode SNP \
	-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 95.0 -tranche 90.0 \
	--rscript-file snps_plots.R \
	--tranches-file snps.tranches \
	-O snps.recal \
	1>log.snprecalibrator 2>&1

gatk --java-options "-Xmx40G -Djava.io.tmpdir=./" ApplyVQSR \
	-R ${ref} \
	-V combined.raw.vcf.gz \
	--ts-filter-level 99.0 \
	--tranches-file snps.tranches \
	--recal-file snps.recal \
	-mode SNP \
	-O snps_VQSR.vcf.gz \
	1>log.applyvqsrsnp 2>&1

gatk --java-options "-Xmx40G -Djava.io.tmpdir=./" VariantRecalibrator \
        -R ${ref} \
        -V snps_VQSR.vcf.gz \
        -resource:mills,known=true,training=true,truth=true,prior=12.0 ${gatkresource}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
        -an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
        -mode INDEL \
	--max-gaussians 6 \
        --rscript-file snps_indels_plots.R \
        --tranches-file snps_indels.tranches \
        -O snps_indels.recal \
        1>log.indelrecalibrator 2>&1

gatk --java-options "-Xmx40G -Djava.io.tmpdir=./" ApplyVQSR \
        -R ${ref} \
        -V snps_VQSR.vcf.gz \
        --ts-filter-level 99.0 \
        --tranches-file snps_indels.tranches \
        --recal-file snps_indels.recal \
        -mode INDEL \
        -O VQSR.vcf.gz \
        1>log.applyvqsrindel 2>&1
```
```sh
nohup bash gatk.sh &
```

### Summary and Visualize .vcf file
- creat index for .vcf files
```sh
$vcf=.vcf file
gatk IndexFeatureFile -I $vcf
```
- summary and visualize
```sh
## VariantQC required, "DISCVRSeq-1.3.8.jar"
gatkresource=path_to_gatk_bundle
java -jar ~/bin/DISCVRSeq-1.3.8.jar  VariantQC -O raw.vcf.html  -R ${gatkresource}/Homo_sapiens_assembly38.fasta -V combined.raw.vcf.gz
java -jar ~/bin/DISCVRSeq-1.3.8.jar  VariantQC -O VQSR.vcf.html  -R ${gatkresource}/Homo_sapiens_assembly38.fasta -V VQSR.vcf.gz
```
