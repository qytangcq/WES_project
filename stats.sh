# ref=~/qytang/gatk/hg38/Homo_sapiens_assembly38.fasta
ref=/mnt/sdafir/qytang/gatkresource/Homo_sapiens_assembly38.fasta
cat 0 | while read id
do
bam=./${id}.bam
samtools stats -@ 16 --reference ${ref} ${bam} > ./stats/${id}.stat
plot-bamstats -p ./stats/${id}  ./stats/${id}.stat 
done
