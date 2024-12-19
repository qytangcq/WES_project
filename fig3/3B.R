library(tidyverse)
#BiocManager::install("maftools")

library(maftools)
library(RColorBrewer) #display.brewer.all() brewer.pal.info
library(ggsci) #pal_npg("nrc", alpha = 0.7)(9)

var.annovar.maf = annovarToMaf(annovar = "data.3BF.txt", 
                               Center = 'NA', 
                               refBuild = 'hg38', 
                               tsbCol = 'Sample', 
                               table = 'refGene',
                               sep = "\t")
write.table(var.annovar.maf,file="all_annovar.maf",quote= F,sep="\t",row.names=F)
pdata <- read.csv("clinical_phenotype") %>%  mutate(outcome=if_else(outcome=="1","Non_bleeding","Bleeding")) %>% rename(Tumor_Sample_Barcode=sample)
var_maf = read.maf(maf ="all_annovar.maf",clinicalData = pdata)


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