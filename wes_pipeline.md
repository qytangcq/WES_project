```sh
mkdir wes_project
wd=path_to_work_directory/wes_project
cd $wd
mkdir {raw,clean,qc,align,mutation,bamqc,annotation} ##with all fastq file saved in raw fold
```

### Quality Control
- summary stats of raw fastq data
```sh
conda activate wes
cd qc && mkdir {raw,clean}
fastqc $wd/raw/* -t 10 -o $wd/qc/raw   ##10 threads
multiqc -n raw $wd/qc/raw              ##summary reports
```

- trim low quality reads, adapters, etc.
```sh
#!/bin/bash
wd=path_to_work_directory/wes_project
source ~/anaconda3/bin/activate wes
for ID in {1..42}
do

trim_galore --paired -q 20 --phred33 --length 36 -e 0.1 --gzip -j 4 -o $wd/clean $wd/raw/${ID}_*R1* $wd/raw/${ID}_*R2*

done
```

```sh
nohup bash qc.sh &
```

- summary stats of clean fastq data
```sh
conda activate wes
fastqc $wd/clean/* -t 10 -o $wd/qc/clean   ##10 threads
multiqc -n clean $wd/qc/clean              ##summary reports
```

### Align to Reference of Human Genome (hg38)

- buid index for human genome reference
```sh
##reference data required by GATK should be prepared
gatkresource=path_to_gatk_bundle
cd ${gatkresource}
conda activate wes
gzip -d Homo_sapiens_assembly38.fasta.gz
bwa index Homo_sapiens_assembly38.fasta ##creating index file .amb .ann .bwt .pac .sa
```
- Do alignment
```sh
#!/bin/bash
gatkresource=path_to_gatk_bundle
INDEX=${gatkresource}/Homo_sapiens_assembly38.fasta
wd=path_to_work_directory/wes_project
source ~/anaconda3/bin/activate wes
cd $wd/align
ls ../clean/*1.fq.gz | sort > 1
ls ../clean/*2.fq.gz | sort > 2
cut -d"/" -f 3 1 | cut -d"_" -f 1 > 0
paste 0 1 2 > config
cat config | while read id
do
arr=($id)
fq1=${arr[1]}
fq2=${arr[2]}
sample=${arr[0]}
bwa mem -t 5 -R "@RG\tID:${sample}\tSM:${sample}\tLB:WGS\tPL:ILLUMINA" $INDEX $fq1 $fq2 | samtools sort -@ 5 -o ${sample}.bam -
done
```
```
nohup bash align.sh &
```
- summary stats of .bam files
```sh
#!/bin/bash
wd=path_to_work_directory/wes_project
gatkresource=path_to_gatk_bundle
INDEX=${gatkresource}/Homo_sapiens_assembly38.fasta
cd $wd/align
source ~/anaconda3/bin/activate wes
cat 0 | while read id
do
bam=./${id}.bam
samtools stats -@ 16 --reference ${INDEX} ${bam} > $wd/bamqc/${id}.stat
plot-bamstats -p $wd/bamqc/${id}  $wd/bamqc/${id}.stat
done
```
```
nohup bash bamqc.sh &
```

### GATK Best Practice to Call Germline Variants
```sh
#!/bin/bash
gatkresource=path_to_gatk_bundle
ref=${gatkresource}/Homo_sapiens_assembly38.fasta
snp=${gatkresource}/dbsnp_146.hg38.vcf.gz
indel=${gatkresource}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
wd=path_to_work_directory/wes_project
source ~/anaconda3/bin/activate wes

cd $wd/mutation

find $wd/align -name "*.bam" | while read id; do basename ${id} .bam; done > config

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
```
```sh
nohup bash gatk.sh &
```

### Summary and Visualize .vcf file
- creat index for .vcf files
```sh
$vcf=.vcf file
gatk IndexFeatureFile -I $vcf
```
- summary and visualize
```sh
## VariantQC required, "DISCVRSeq-1.3.8.jar"
gatkresource=path_to_gatk_bundle
java -jar ~/bin/DISCVRSeq-1.3.8.jar  VariantQC -O raw.vcf.html  -R ${gatkresource}/Homo_sapiens_assembly38.fasta -V combined.raw.vcf.gz
java -jar ~/bin/DISCVRSeq-1.3.8.jar  VariantQC -O VQSR.vcf.html  -R ${gatkresource}/Homo_sapiens_assembly38.fasta -V VQSR.vcf.gz
```

Individuals' WES data from both case and control group should be processed according to the same protocol, joint-called together or in a group.
For example, **case.VQSR.vcf.gz** and **control.VQSR.vcf.gz**
Then, downstream analysis should be undertaken.

### Downstream Analysis of .vcf files
take **case.VQSR.vcf.gz** as example

- normalization of variants
```sh
gatkresource=path_to_gatk_bundle
conda activate wes
bcftools norm -m-both case_VQSR.vcf.gz -o  case_tmp.vcf
bcftools norm -f ${gatkresource}/Homo_sapiens_assembly38.fasta case_tmp.vcf -o case_VQSR_norm.vcf
```
- filtration of variants, by missing rate/minor allele count/minimum phred-sclaed Quality
```sh
conda activate wes
vcftools --vcf case_VQSR_norm.vcf --max-missing 0.5 --mac 3 --minQ 30  --remove-filtered VQSRTrancheSNP99.90to100.00  --remove-filtered VQSRTrancheSNP99.90to100.00+ --recode --recode-INFO-all --out case_VQSR_norm_filtered
mv case_VQSR_norm_filtered.recode.vcf case_VQSR_norm_filtered.vcf
```

- (optional) annotate variants with rsid (dbSNP database)
```sh
##SnpSift required, should be installed according to official guidance, so as its ref data (00-All.vcf.gz)
SnpSift=path_to_SnpSift.jar
ref=path_to_00-All.vcf.gz
java -jar ${SnpSift} annotate -id ${ref} case_VQSR_norm_filtered.vcf > case_VQSR_norm_filtered_snpID.vcf
```
- (optional) annotate variants with chr:location as SNP's name (for following analysis in case of that rsid could not be applicable for all variants detected)
```
cat case_VQSR_norm_filtered.vcf |awk -v OFS="\t" '/^#/|| $3=$1"_"$2{print $0}' > case_VQSR_norm_filtered_snpID.vcf
```

#### SNV-based association test between phenotype and variants, using variants' allele frequnecy
- annotate .vcf with ANNOVAR
```
###with ANNOVAR installed and its ref data downloaded, directory written to $PATH
wd=path_to_work_directory/wes_project
db=path_to_annovar_humandb
convert2annovar.pl -format vcf4 -allsample -withfreq case_VQSR_norm_filtered_snpID.vcf  -out $wd/annotation/case.avinput
table_annovar.pl $wd/annotation/case.avinput ${db} -buildver hg38 -out $wd/annotation/case -remove -protocol refGene,avsnp150,clinvar_20210501,dbnsfp42c,exac03,intervar_20180118,EAS.sites.2015_08 -operation g,f,f,f,f,f,f -nastring NA -csvout
```

- annotation data process, for case and control group separately
```
##requiring csvtk
##case group
cd $wd/annotation
conda activate wes

cat case.hg38_multianno.csv | csvtk filter2 -f '$9!="synonymous SNV" && $9!="NA" && $9!="unknown"' | csvtk filter2 -f '$ExAC_EAS<0.01' >case.unsyno.1%.csv #filter out synonymous variants and keep rare Variants

cat case.hg38_multianno.csv | csvtk filter2 -f '$9 == "synonymous SNV" || $9 == "NA" || $9 == "unknown"' | csvtk filter2 -f '$144<0.01' | csvtk filter2 -f '$SIFT_pred == "D"' > case.syno.1%.D.csv #rescue the rare synonymous variants predicted to be deleterious by SIFT algorithm

csvtk concat case.syno.1%.D.csv  case.unsyno.1%.csv > case.D.csv

cat annotate/case.avinput |tr " \t" "," >case.avinput.csv

##control group
cat control.hg38_multianno.csv | csvtk filter2 -f '$9!="synonymous SNV"' | csvtk filter2 -f '$ExAC_EAS<0.01' >control.unsyno.1%.csv  #filter out synonymous variants and keep rare Variants

cat control.hg38_multianno.csv | csvtk filter2 -f '$9 == "synonymous SNV" || $9 == "NA" || $9 == "unknown"' | csvtk filter2 -f '$144<0.01' | csvtk filter2 -f '$SIFT_pred == "D"' > control.syno.1%.D.csv #rescue the rare synonymous variants predicted to be deleterious by SIFT algorithm

csvtk concat control.syno.1%.D.csv  control.unsyno.1%.csv > control.D.csv

cat annotate/control.avinput |tr " \t" "," >control.avinput.csv

##intersection of cases' and controls' variants
csvtk inter -f 1,2,3,4,5 bleeding.D.csv  control.D.csv > case_control.harmonised.csv

##add allele frequency information in case and control group
csvtk -H join -f 1,2,3,4,5 case_control.harmonised.csv case.avinput.csv | csvtk -H cut -f 1,2,3,4,5,6 |csvtk -H join -f 1,2,3,4,5 - control.avinput.csv |csvtk -H cut -f 1,2,3,4,5,6,7 | csvtk -H add-header -n chr,start,end,ref,alt,af_case,af_control > case_control.harmonised.af.csv
```

```R
##wd=path_to_work_directory/wes_project

rm(list=ls())
setwd("$wd/annotation")
library("tibble")
library(tidyr)

##load data
data <- read.csv("harmonized_anno_af",header = T,sep = ",")
data <- data[,c(1,2,4,5,145,146)]
data$Start <- as.character(data$Start)
data <- unite(data,id,Chr,Start,Ref,Alt,sep=":")

data$ac_case <- round(data$af_case*45*2,digits = 0)
data$ac_control <- round(data$af_control*40*2,digits = 0)
data <- data[,c(1,4,5)]

data1 <- data
data1$ac_case <- 45*2-data1$ac_case
data1$ac_control <- 40*2-data1$ac_control

data_snp <- rbind(data,data1)

tmp<-data.frame()
for (i in unique(data_snp$id)){
  sub_data<-subset(data_snp,data_snp$id==i)
  mat <- as.matrix(sub_data[,c(2,3)])
  output_fisher_snp<-chisq.test(mat) ##optional:fisher.test
  output_fisher_snp<-data.frame(output_fisher_snp$p.value)
  output_fisher_snp<-add_column(output_fisher_snp,ID=i)
  output_fisher_snp<-add_column(output_fisher_snp,ac_case=sub_data[,2][1])
  output_fisher_snp<-add_column(output_fisher_snp,ac_control=sub_data[,3][1])
  tmp<-rbind(tmp,output_fisher_snp)
}

tmp$p_ad<- p.adjust (tmp$output_fisher_snp.p.value, method="BH")

tmp <- tmp %>% separate(ID,into = c("Chr","Position","Ref","Alt"),sep = ":",remove = T)

tmp <- tmp[,c(2:8)]
write.csv(tmp,file = "fisher_snp_af_detail.csv",quote = F,row.names = F)
```
```sh
conda activate wes
cat fisher_snp_af_detail.csv | csvtk filter2 -f '$p_ad<0.05' |csvtk mutate -f Position -n end |csvtk cut -f 1,2,6,3,4,5|csvtk del-header |csvtk -H add-header -n Chr,Start,End,Ref,Alt,p_adj | csvtk join -f 1,2,3,4,5 - case.hg38_multianno.csv > fisher.sig.snp.annonated.csv
```

#### GWAS (optional in case of low power when sample size is small)
```sh
##requiring plink
wd=path_to_work_directory/wes_project
cd $wd
mkdir gwas && cd gwas
conda activate wes

bcftools merge -m both case_VQSR_norm_filtered.vcf control_VQSR_norm_filtered.vcf -o combined.vcf ##combine case and control vcf

plink --vcf combined.vcf --recode --out var --double-id --allow-extra-chr
plink --file var --make-bed --out var --allow-extra-chr --allow-no-sex ## trans .vcf file into plink binary format

plink --bfile var--missing --allow-extra-chr --allow-no-sex
Rscript --no-save hist_miss.R
plink --bfile var --geno 0.3 --make-bed --out var_2 --allow-extra-chr  --allow-no-sex
plink --bfile var_2 --mind 0.3 --make-bed --out var_3 --allow-extra-chr --allow-no-sex
plink --bfile var_3 --check-sex -allow-extra-chr  --allow-no-sex
Rscript --no-save gender_check.R
plink --bfile var_3 --impute-sex --make-bed --out var_4 -allow-extra-chr
awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' var_4.bim > snp_1_22.txt
plink --bfile var_4 --extract snp_1_22.txt --make-bed --out var_5 --allow-extra-chr  --allow-no-sex
plink --bfile var_5 --freq --out MAF_check --allow-extra-chr  --allow-no-sex
Rscript --no-save MAF_check.R
plink --bfile var_5 --maf 0.05 --make-bed --out var_6 --allow-extra-chr --allow-no-sex
plink --bfile var_6 --hardy --allow-extra-chr --allow-no-sex
awk '{ if ($9 <0.00001) print $0 }' plink.hwe>plinkzoomhwe.hwe
plink --bfile var_6 --hwe 1e-6 --make-bed --out var_hwe_filter_step1 --allow-extra-chr --allow-no-sex
plink --bfile var_hwe_filter_step1 --hwe 1e-10 --hwe-all --make-bed --out var_7  --allow-extra-chr --allow-no-sex
plink --bfile var_7 --exclude inversion.txt --range --indep-pairwise 50 5 0.2 --out indepSNP --allow-extra-chr --allow-no-sex ##inversion.txt: genome region with high variance
plink --bfile var_7 --extract indepSNP.prune.in    --het --out R_check --allow-extra-chr --allow-no-sex
Rscript --no-save check_heterozygosity_rate.R
Rscript --no-save heterozygosity_outliers_list.R
sed 's/"// g' fail-het-qc.txt | awk '{print$1, $2}'> het_fail_ind.txt
plink --bfile var_7 --remove het_fail_ind.txt --make-bed --out var_8 --allow-extra-chr --allow-no-sex
plink --bfile var_8 --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2 --allow-extra-chr --allow-no-sex
awk '{ if ($8 >0.9) print $0 }' pihat_min0.2.genome>zoom_pihat.genome
Rscript --no-save Relatedness.R
plink --bfile var_8 --filter-founders --make-bed --out var_9 --allow-extra-chr --allow-no-sex


plink --bfile var_7 --pca 10 --out var_pca --allow-extra-chr  --allow-no-sex ##pricipal component analysis

#BOFORE association test, .fam file should be edit to add phenotype information

##Association tests
plink --bfile var_7 --assoc  --adjust --out result_gwas_asso_adjust --covar var_pca.eigenvec --allow-extra-chr --allow-no-sex
plink --bfile var_7 --fisher --allow-extra-chr --allow-no-sex --out result_gwas_fisher
plink --bfile var_7 --logistic --covar bleeding_pca.eigenvec --hide-covar --allow-extra-chr --allow-no-sex --out result_gwas_logistic

cat result_gwas_logistic.assoc.logistic |sed 's/[ ][ ]*/,/g' |sed 's/,//' |csvtk -t del-header |csvtk -H add-header -n chr,snp,bp,a1,test,nmiss,or,stat,p |csvtk filter2 -f '$test=="ADD"' > result_gwas_logistic.assoc.logistic.add

cat result_gwas_fisher.assoc.fisher |csvtk cut -f 2,3,4,5,6,7|csvtk filter2 -f '$P<0.05' | csvtk mutate -f BP -n end | csvtk cut -f 1,2,7,3,4,5,6 |csvtk del-header |csvtk -H add-header -n Chr,Start,End,Ref,Alt,OR,p_adj | csvtk join -f 1,2,3,4,5 - bleeding.hg38_multianno.csv > gwas_fisher_sig_snp_annotate.csv
```

#### Gene-based association test between phenotype and variants with SKAT

In our study, previous sreening of variants led to a group of candidate SNV, annotated as **case_control.harmonised.af.csv**

based on case_control.harmonised.af.csv, build .bed files recording the chrom,chromStart and chromEnd information, so as a subset of vcf file could be created
```sh
wd=path_to_work_directory/wes_project
cd $wd && mkdir SKAT && cd SKAT
conda activate wes

cat $wd/annotation/case_control.harmonised.ad.csv |csvtk cut -f 1,2,3,7 > SNP_list

mkdir genebed
```
```R
#wd=path_to_work_directory/wes_project
setwd(dir="$wd/SKAT")

library(tidyverse)

#process of snp file
snp_list <- read.csv("SNP_list",header = T)
gene_name <- unique(snp_list$Gene.refGene)

for (i in gene_name) {
  sub_list <- subset(snp_list,Gene.refGene==i)[,1:3]
  colnames(sub_list) <- c("chrom","chromStart","chromEnd")
  write.table(sub_list,file = paste0("./genebed/",i,".bed"),quote = F,sep = "\t",col.names = T,row.names = F)
}
```
```sh
cd $wd/SKAT
cat $wd/annotation/case_control.harmonised.ad.csv |csvtk cut -f 1,2,3 |csvtk del-header |csvtk -H add-header -n chrom,chromStart,chromEnd |tr "," "\t" >all.bed

vcftools --vcf $wd/gwas/combined.vcf --bed all.bed --recode  --out selected1322snp

cat SNP_list | csvtk cut -f 4 |sed '1d' >gene_name

mkdir genevcf
cat gene_name |while read i;do (vcftools --vcf selected1322snp.recode.vcf  --bed genebed/${i}.bed --recode --out genevcf/${i} );done

mkdir plink
cat gene_name |while read i;do (plink --vcf genevcf/${i}.recode.vcf recode --out plink/${i} --double-id --allow-extra-chr);done
cat gene_name |while read i;do (plink --file plink/${i} --make-bed --out plink/${i} --allow-extra-chr --allow-no-sex);done

##with plink2
cat gene_name |while read i;do (plink2 --bfile plink/${i} --make-bed out plink/${i}_undup --rm-dup force-first --allow-extra-chr --allow-no-sex);done
```
```R
#wd=path_to_work_directory/wes_project
setwd(dir="$wd/SKAT")

library(tidyverse)
library(SKAT)

snp_list <- read.csv("SNP_list",header = T)
gene_name <- unique(snp_list$Gene.refGene)

re <- data.frame()

for (i in gene_name) {
  i="FCGBP"
  file.bed <- paste0("./plink/",i,"_undup.bed")
  file.bim <- paste0("./plink/",i,"_undup.bim")
  file.fam <- paste0("./plink/",i,"_undup.fam")
  file.SSD <- paste0("./plink/",i,"_undup.SSD")
  file.info <- paste0("./plink/",i,"_undup.info")

  bim <- read.table(file.bim)[,1:2] %>% mutate(v3=i)
  write.table(bim[,c(3,2)],file = paste0("./plink/",i,"_undup.SetID"),col.names = F,row.names = F,quote = F,sep = "\t")
  file.setID <- paste0("./plink/",i,"_undup.SetID")
  Generate_SSD_SetID(file.bed,file.bim,file.fam,file.setID,file.SSD,file.info)
  FAM <- Read_Plink_FAM(file.fam,Is.binary = F)
  y <- FAM$Phenotype-1

  SSD.INFO <- Open_SSD(file.SSD,file.info)
  SSD.INFO$nSample
  SSD.INFO$nSets

  obj <- SKAT_Null_Model(y~1,out_type = "D")

  out <- SKAT_CommonRare.SSD.All(SSD.INFO,obj,method="C")
  Close_SSD()

  p <- out$results$P.value
  tmp <- data.frame(gene=i,p_value=p)
  re <- rbind(re,tmp)
}

write.csv(re,file="./results_SKAT_Comonrare.csv")
```

##### (optional) waterfall plot, lollipop plot etc

```sh
bcftools query -l selected1322snp_snpID.vcf > indiv.list
mkdir indivvcf
cat indiv.list |while read i;do (touch indivvcf/${i}.list;echo ${i}>> indivvcf/${i}.list);done
cat indiv.list |while read i;do (bcftools view -S indivvcf/${i}.list selected1322snp_snpID.vcf >indivvcf/${i}.vcf);done

cat indiv.list  |while read i;do (cat indivvcf/${i}.vcf |grep -v "0/0" |grep -v "0|0" >indivvcf/${i}.recode.vcf);done

db=path_to_annovar_humandb
mkdir indivanno
cat indiv.list |while read i;do (table_annovar.pl indivvcf/${i}.recode.vcf ${db} -buildver hg38 -out ./indivanno/${i} -remove -protocol refGene -operation g -nastring . -vcfinput -polish);done

cd indivanno  && for i in *.hg38_multianno.txt;do (sample=`echo $i|awk -F '.' '{print $1}'`;cut -f '1-10' ${i}|sed '1d'|sed "s/$/\t${sample}/">>$wd/SKAT/all_anno.txt);done
```
```R
#wd=path_to_work_directory/wes_project
setwd(dir="$wd/SKAT")
library(tidyverse)
library(maftools)
library(RColorBrewer) #display.brewer.all() brewer.pal.info
library(ggsci) #pal_npg("nrc", alpha = 0.7)(9)

var.annovar.maf = annovarToMaf(annovar = "./all_anno.txt",
                               Center = 'NA',
                               refBuild = 'hg38',
                               tsbCol = 'Sample',
                               table = 'refGene',
                               sep = "\t")
write.table(var.annovar.maf,file="all_annovar.maf",quote= F,sep="\t",row.names=F)

##clinical data required to set group
pdata <- read.csv("clinical_phenotype") %>%  mutate(outcome=if_else(outcome=="1","Non_bleeding","Bleeding")) %>% rename(Tumor_Sample_Barcode=sample)
var_maf = read.maf(maf ="all_annovar.maf",clinicalData = pdata)

pdf("summary.pdf", width=6, height=6)
plotmafSummary(maf = var_maf, rmOutlier = TRUE, addStat = 'median',)
dev.off()

gene <- read.csv("results_SKAT_Comonrare.csv")  %>%  as_tibble() %>% subset(p_value<0.01)

col_outcome=pal_jama("default", alpha = 0.8)(2)
assign_outcome=setNames(col_outcome,unique(pdata$outcome))
colors <- pal_npg("nrc", alpha = 0.8)(9)
names(colors) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del',
  'Nonstop_Mutation'
)

pdf("All.oncoplot.clin.pdf", width=8, height=6)
oncoplot(maf =var_maf, fontSize = 0.45,showTumorSampleBarcodes = F,
         SampleNamefontSize=0.7,titleFontSize=1.2,
         legendFontSize=1,removeNonMutated=T,
         writeMatrix=T,draw_titv = T,sortByMutation = T,
         genes = gene$gene,keepGeneOrder = F,
         clinicalFeatures = "outcome",sortByAnnotation = T, anno_height = 0.5,
         annotationColor = list(outcome=assign_outcome),
         colors = colors,bgCol = "white") #
dev.off()

var_maf_backup=var_maf
var_maf@data[which(Hugo_Symbol=="FCGBP" & is.na(Variant_Type) & Variant_Classification=="In_Frame_Ins"),"Variant_Type"]="INS"
pdf("lollipopplot_FCGBP.pdf", width=8, height=6)
lollipopPlot(maf = var_maf, gene = 'FCGBP', AACol = "aaChange",
             showMutationRate = F,labelPos = "1221",repel = F,showDomainLabel = F,printCount = T)
dev.off()

install.packages("wordcloud")
library(wordcloud)
devtools::source_gist(id = "https://gist.github.com/PoisonAlien/3f8752a89c1d63f64afe55b441621223")
pdf("genecloud.pdf", width=20, height=20)
geneCloud(input = var_maf_backup, top=200)
dev.off()
```
