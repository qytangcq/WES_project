# exon_bed=~/qytang/gatk/hg38/hg38.exon.bed
exon_bed=/mnt/sdafir/qytang/gatkresource/hg38.exon.bed
ls *_bqsr.bam | while read id
do
qualimap bamqc --java-mem-size=20G -gff $exon_bed -bam $id
done
