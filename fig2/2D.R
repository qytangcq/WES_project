library(tidyverse)
#BiocManager::install("maftools")

library(maftools)
library(RColorBrewer) #display.brewer.all() brewer.pal.info
library(ggsci) #pal_npg("nrc", alpha = 0.7)(9)

var.annovar.maf = annovarToMaf(annovar = "dat.2D.txt", 
                               Center = 'NA', 
                               refBuild = 'hg38', 
                               tsbCol = 'Sample', 
                               table = 'refGene',
                               sep = "\t")
write.table(var.annovar.maf,file="all_annovar.maf",quote= F,sep="\t",row.names=F)
pdata <- read.csv("clinical_phenotype") %>%  mutate(outcome=if_else(outcome=="1","Non_bleeding","Bleeding")) %>% rename(Tumor_Sample_Barcode=sample)
var_maf = read.maf(maf ="all_annovar.maf",clinicalData = pdata)


pdf("summary.pdf", width=6, height=6)
plotmafSummary(maf = var_maf, rmOutlier = TRUE, addStat = 'median',dashboard = T)
dev.off()
