INDEX=/mnt/sdafir/qytang/gatkresource/Homo_sapiens_assembly38.fasta
export SENTIEON_LICENSE=/mnt/sdafir/qytang/sentieon/Southeast_University_eval.lic

for i in {104..105}
do
sentieon driver -t 40 -r ${INDEX} \
  -i /mnt/sdafir/qytang/bleeding_project/dedup_bam/${i}_dedup.bam --algo QualCal ${i}_RECAL_DATA.TABLE

sentieon driver -t 40 -r ${INDEX} -i /mnt/sdafir/qytang/bleeding_project/dedup_bam/${i}_dedup.bam \
  -q ${i}_RECAL_DATA.TABLE --algo QualCal ${i}_RECAL_DATA.TABLE.POST \
  --algo ReadWriter ${i}_dedup_recalibrated.bam
sentieon driver -t 40 --algo QualCal --plot \
  --before ${i}_RECAL_DATA.TABLE --after ${i}_RECAL_DATA.TABLE.POST ${i}_RECAL_RESULT.CSV
sentieon plot QualCal -o ${i}_BQSR_PDF ${i}_RECAL_RESULT.CSV

done

for i in {108..145}
do
sentieon driver -t 40 -r ${INDEX} \
  -i /mnt/sdafir/qytang/bleeding_project/dedup_bam/${i}_dedup.bam --algo QualCal ${i}_RECAL_DATA.TABLE

sentieon driver -t 40 -r ${INDEX} -i /mnt/sdafir/qytang/bleeding_project/dedup_bam/${i}_dedup.bam \
  -q ${i}_RECAL_DATA.TABLE --algo QualCal ${i}_RECAL_DATA.TABLE.POST \
  --algo ReadWriter ${i}_dedup_recalibrated.bam
sentieon driver -t 40 --algo QualCal --plot \
  --before ${i}_RECAL_DATA.TABLE --after ${i}_RECAL_DATA.TABLE.POST ${i}_RECAL_RESULT.CSV
sentieon plot QualCal -o ${i}_BQSR_PDF ${i}_RECAL_RESULT.CSV

done

sentieon driver -t 40 -r ${INDEX} \
         -i 104_dedup_recalibrated.bam  -i 105_dedup_recalibrated.bam \
         $(for y in {108..145};do echo "-i ${y}_dedup_recalibrated.bam";done) \
        --algo Genotyper sentieon_control.vcf.gz

gatkresource=/mnt/sdafir/qytang/gatkresource

sentieon driver -t 40  -r ${INDEX} --algo VarCal -v sentieon_control.vcf.gz \
        --annotation DP --annotation QD --annotation FS --annotation SOR \
        --annotation ReadPosRankSum --annotation MQRankSum \
        --tranches_file snps_control.tranches \
        --resource ${gatkresource}/hapmap_3.3.hg38.vcf.gz --resource_param hapmap,known=false,training=true,truth=true,prior=15.0 \
        --resource ${gatkresource}/1000G_omni2.5.hg38.vcf.gz --resource_param omini,known=false,training=true,truth=false,prior=12.0 \
        --resource ${gatkresource}/1000G_phase1.snps.high_confidence.hg38.vcf.gz --resource_param 1000G,known=false,training=true,truth=false,prior=10.0 \
        --resource ${gatkresource}/dbsnp_146.hg38.vcf.gz --resource_param dbsnp,known=true,training=false,truth=false,prior=2.0 \
        --var_type SNP  \
        --tranche 100.0 --tranche 99.9 --tranche 99.0 --tranche 95.0 --tranche 90.0 \
        snps_control.recal \
        1>log.snprecalibrator. 2>&1


sentieon driver -t 40  -r ${INDEX} --algo VarCal -v sentieon_control.vcf.gz \
        --resource ${gatkresource}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --resource_param mills,known=true,training=true,truth=true,prior=12.0\
        --annotation DP --annotation QD --annotation FS --annotation SOR \
        --annotation ReadPosRankSum --annotation MQRankSum \
        --tranches_file indels_control.tranches \
        --var_type INDEL \
        --max_gaussians 6 \
        indels_control.recal \
        1>log.indelrecalibrator 2>&1

sentieon driver -t 40   -r ${INDEX} --algo ApplyVarCal -v sentieon_control.vcf.gz \
        --recal snps_control.recal --tranches_file snps_control.tranches \
        --var_type SNP --sensitivity 99.0 \
        sentieon_control_snps_VQSR.vcf.gz \
        1>log.applyvqsrsnp 2>&1

sentieon driver -t 40   -r ${INDEX} --algo ApplyVarCal -v  sentieon_control.vcf.gz \
        --recal indels_control.recal --tranches_file indels_control.tranches \
        --var_type INDEL --sensitivity 99.0 \
        sentieon_control_indel_VQSR.vcf.gz \
        1>log.applyvqsrindel 2>&1




