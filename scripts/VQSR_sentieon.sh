gatkresource=/mnt/sdafir/qytang/gatkresource
INDEX=/mnt/sdafir/qytang/gatkresource/Homo_sapiens_assembly38.fasta
export SENTIEON_LICENSE=/mnt/sdafir/qytang/sentieon/Southeast_University_eval.lic


sentieon driver -t 40  -r ${INDEX} --algo VarCal -v sentieon_combined.vcf.gz \
	--annotation DP --annotation QD --annotation FS --annotation SOR \
	--annotation ReadPosRankSum --annotation MQRankSum \
	--tranches_file snps.tranches \
	--resource ${gatkresource}/hapmap_3.3.hg38.vcf.gz --resource_param hapmap,known=false,training=true,truth=true,prior=15.0 \
	--resource ${gatkresource}/1000G_omni2.5.hg38.vcf.gz --resource_param omini,known=false,training=true,truth=false,prior=12.0 \
        --resource ${gatkresource}/1000G_phase1.snps.high_confidence.hg38.vcf.gz --resource_param 1000G,known=false,training=true,truth=false,prior=10.0 \
        --resource ${gatkresource}/dbsnp_146.hg38.vcf.gz --resource_param dbsnp,known=true,training=false,truth=false,prior=2.0 \
	--var_type SNP	\
	--tranche 100.0 --tranche 99.9 --tranche 99.0 --tranche 95.0 --tranche 90.0 \
	snps.recal \
	1>log.snprecalibrator 2>&1


sentieon driver -t 40  -r ${INDEX} --algo VarCal -v sentieon_combined.vcf.gz \
        --resource ${gatkresource}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --resource_param mills,known=true,training=true,truth=true,prior=12.0\
	--annotation DP --annotation QD --annotation FS --annotation SOR \
        --annotation ReadPosRankSum --annotation MQRankSum \
	--tranches_file indels.tranches \
        --var_type INDEL \
        --max_gaussians 6 \
        indels.recal \
        1>log.indelrecalibrator 2>&1

sentieon driver -t 40   -r ${INDEX} --algo ApplyVarCal -v sentieon_combined.vcf.gz \
	--recal snps.recal --tranches_file snps.tranches \
	--var_type SNP --sensitivity 99.0 \
	snps_VQSR.vcf.gz \
	1>log.applyvqsrsnp 2>&1 

sentieon driver -t 40   -r ${INDEX} --algo ApplyVarCal -v snps_VQSR.vcf.gz \
        --recal indels.recal --tranches_file indels.tranches \
        --var_type INDEL --sensitivity 99.0 \
        VQSR_sentieon.vcf.gz \
	1>log.applyvqsrindel 2>&1









