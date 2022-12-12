# INDEX=/home/qytang/qytang/gatk/hg38/bwa_index/Homo_sapiens_assembly38.fasta
INDEX=/mnt/sdafir/qytang/gatkresource/Homo_sapiens_assembly38.fasta
export SENTIEON_LICENSE=/mnt/sdafir/qytang/sentieon/Southeast_University_eval.lic

ls ../trim/*1.fq.gz | sort > 1
ls ../trim/*2.fq.gz | sort > 2
cut -d"/" -f 3 1 | cut -d"_" -f 1 > 0
paste 0 1 2 > config
cat config | while read id
do
arr=($id)
fq1=${arr[1]}
fq2=${arr[2]}
sample=${arr[0]}

(sentieon bwa mem -R "@RG\tID:$sample\tSM:$sample\tLB:WES\tPL:ILLUMINA" \
	-t 40 $INDEX $fq1 $fq2 || echo -n 'error') \
	| sentieon util sort -r $INDEX -o $sample.bam -t 30 --sam2bam -i -
#bwa mem -t 5 -R "@RG\tID:$sample\tSM:$sample\tLB:WGS\tPL:ILLUMINA" $INDEX $fq1 $fq2 | samtools sort -@ 5 -o $sample.bam -
done
