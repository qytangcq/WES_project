# ref=~/qytang/gatk/hg38/Homo_sapiens_assembly38.fasta
# snp=~/qytang/gatk/hg38/dbsnp_146.hg38.vcf.gz
# indel=~/qytang/gatk/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

ref=/mnt/sdafir/qytang/gatkresource/Homo_sapiens_assembly38.fasta
snp=/mnt/sdafir/qytang/gatkresource/dbsnp_146.hg38.vcf.gz
indel=/mnt/sdafir/qytang/gatkresource/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

find -name "*.bam" | while read id; do basename ${id} .bam; done > config

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

#gatk --java-options "-Xmx40G -Djava.io.tmpdir=./" HaplotypeCaller \
#       -R ${ref} \
#       -I ${sample}_bqsr.bam \
#        --dbsnp ${snp} \
#        -O ${sample}_raw.vcf \
#        1>${sample}.log.haplotypecaller 2>&1

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
gatkresource=/mnt/sdafir/qytang/gatkresource

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





	  
