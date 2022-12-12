export OMP_NUM_THREADS=10

excl=/mnt/sdafir/qytang/bleeding_project/sv/human.hg38.excl.tsv
genome=/mnt/sdafir/qytang/gatkresource/Homo_sapiens_assembly38.fasta

delly call -t ALL -g $genome -x $excl -o 1_sv.bcf /mnt/sdafir/qytang/bleeding_project/align/1.bam


