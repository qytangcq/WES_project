##waterfall plot
##TANG QIuyi
##20221205


library(tidyverse)
#BiocManager::install("maftools")

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
pdata <- read.csv("clinical_phenotype") %>%  mutate(outcome=if_else(outcome=="1","Non_bleeding","Bleeding")) %>% rename(Tumor_Sample_Barcode=sample)
var_maf = read.maf(maf ="all_annovar.maf",clinicalData = pdata)



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

var_maf_backup=var_maf
var_maf@data[which(Hugo_Symbol=="FCGBP" & is.na(Variant_Type) & Variant_Classification=="In_Frame_Ins"),"Variant_Type"]="INS"


###creat waterfall plot with ComplexHeatmap
#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(reshape2)
rm(list = ls())
matMut <- read.table("onco_matrix.txt", header = T, check.names = F, sep = "\t")
matMuttmp = matMut
matMuttmp$gene = row.names(matMuttmp)
mat_long <- melt(matMuttmp, id.vars = "gene", value.name = "Variant_Classification")
levels(factor(mat_long$Variant_Classification))

alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(3, "mm"), h-unit(1, "mm"), 
              gp = gpar(fill = "white", col = NA))
  },
  Frame_Shift_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(3, "mm"), h-unit(1, "mm"), 
              gp = gpar(fill = colors["Frame_Shift_Del"], col = NA))
  },
  Frame_Shift_Ins = function(x, y, w, h) {
    grid.rect(x, y, w-unit(3, "mm"), h-unit(1, "mm"),  
              gp = gpar(fill = colors["Frame_Shift_Ins"], col = NA))
  },
  In_Frame_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(3, "mm"), h-unit(1, "mm"), 
              gp = gpar(fill = colors["In_Frame_Del"], col = NA))
  },
  Missense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(3, "mm"), h-unit(1, "mm"), 
              gp = gpar(fill = colors["Missense_Mutation"], col = NA))
  },
  Multi_Hit = function(x, y, w, h) {
     grid.rect(x, y, w-unit(3, "mm"),h-unit(1, "mm"),
              gp = gpar(fill = colors["Multi_Hit"], col = NA))
  }
  #,
  #Nonstop_Mutation = function(x, y, w, h) {
  #  grid.rect(x, y, w-unit(3, "mm"), h-unit(0.5, "mm"), 
  #            gp = gpar(fill = col["Nonstop_Mutation"], col = NA))
  #},
  #In_Frame_Ins = function(x, y, w, h) {
  #  grid.rect(x, y, w-unit(3, "mm"), h-unit(0.5, "mm"), 
  #            gp = gpar(fill = col["In_Frame_Ins"], col = NA))
  #},
  #Nonsense_Mutation = function(x, y, w, h) {
  #  grid.rect(x, y, w-unit(3, "mm"), h-unit(0.5, "mm"), 
  #            gp = gpar(fill = col["Nonsense_Mutation"], col = NA))
  #}
)

heatmap_legend_param <- list(title = "Alternations", 
                             at = c("Frame_Shift_Del", "Frame_Shift_Ins","In_Frame_Del", "Missense_Mutation","Multi_Hit","Nonstop_Mutation","In_Frame_Ins","Nonsense_Mutation"), 
                             labels = c("Frame_Shift_Del", "Frame_Shift_Ins","In_Frame_Del", "Missense_Mutation","Multi_Hit","Nonstop_Mutation","In_Frame_Ins","Nonsense_Mutation"))


#colors defined before
#col <- c(Frame_Shift_Del = "purple", Frame_Shift_Ins = "blue", Multi_Hit = "orange", In_Frame_Del = "black",Missense_Mutation="green",Nonstop_Mutation="gray",In_Frame_Ins="white",Nonsense_Mutation="yellow")

pdata <- read.csv("clinical_phenotype")
ha <- HeatmapAnnotation(phenotype = as.factor(pdata$outcome), show_annotation_name = TRUE,
                        annotation_name_gp = gpar(fontsize = 7)) ##clinical data

column_title <- ""

oncoPrint(matMut, alter_fun = alter_fun, col = colors, alter_fun_is_vectorized = FALSE)

oncoplot_anno<- oncoPrint(matMut,
          bottom_annotation = ha, #注释信息在底部
          #   top_annotation=top_annotation,
          #right_annotation=NULL,
          alter_fun = alter_fun, 
          col = colors,  
          column_title = "", 
          heatmap_legend_param = heatmap_legend_param,
          row_names_side = "left",
          pct_side = "right",
          #column_order=sample_order,
          #       column_split=3
          alter_fun_is_vectorized = FALSE
)
draw(oncoplot_anno, annotation_legend_side = "left" )
