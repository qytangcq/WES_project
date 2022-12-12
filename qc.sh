
dir=/mnt/sdafir/qytang/bleeding_project/
for i in {1..42}
do
trim_galore --paired -q 20 --phred33 --length 36 -e 0.1 --gzip -j 4 -o $dir ../raw/${i}_*R1* ../raw/${i}_*R2* 
done
