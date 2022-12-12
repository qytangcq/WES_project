ref=/mnt/sdafir/qytang/gatkresource/Homo_sapiens_assembly38.fasta

find ./gvcf_gatk -name "*.g.vcf.gz" > gvcf.list
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











	

