INDEX=/mnt/sdafir/qytang/gatkresource/Homo_sapiens_assembly38.fasta
export SENTIEON_LICENSE=/mnt/sdafir/qytang/sentieon/Southeast_University_eval.lic

for i in {1..48}
do
sentieon driver -t 40 -i ../align/${i}.bam \
  --algo LocusCollector --fun score_info ${i}_SCORE.gz

sentieon driver -t 40 -i ../align/${i}.bam \
  --algo Dedup --rmdup --score_info ${i}_SCORE.gz \
  --metrics ${i}_DEDUP_METRIC_TXT ${i}_dedup.bam
done

for i in {101..145}
do
sentieon driver -t 40 -i ../align/${i}.bam \
  --algo LocusCollector --fun score_info ${i}_SCORE.gz

sentieon driver -t 40 -i ../align/${i}.bam \
  --algo Dedup --rmdup --score_info ${i}_SCORE.gz \
  --metrics ${i}_DEDUP_METRIC_TXT ${i}_dedup.bam
done
