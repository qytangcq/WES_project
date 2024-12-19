
library(RNAseqStat2)
library(tidyverse)
library(ggthemes)
library(EnvStats)
library(ggpubr)
library(patchwork)
library(ClusterGVis)
library(org.Hs.eg.db)
library(scCustomize)

##RNAseq_analysis
counts <- read.csv("data.5BC.csv",header = TRUE)
counts <- counts  %>% column_to_rownames("gene_name")
group <- c(rep("Control",3),rep("Si",3))
data_i <- Create_DEGContainer(species = "Human",
                              dataType = "Counts",
                              idType = "SYMBOL",
                              expMatrix = counts,
                              groupInfo = group,
                              caseGroup = "Si")
data_o <- runALL(object = data_i,dir = "results/DEG",GO = FALSE, KEGG = FALSE)


#### plot DEG and enrichment analysis
TPM <- openxlsx::read.xlsx("TPM_table.xlsx",sheet = "gene")
TPM <- TPM %>% dplyr::filter(gene_type %in% "protein_coding") %>% dplyr::select(gene_name,S4_1,S4_2,S4_3,NC_1,NC_2,NC_3) %>% distinct(.,gene_name,.keep_all = TRUE) %>% column_to_rownames("gene_name")
colnames(TPM) <- c("Si_1","Si_2","Si_3","NC_1","NC_2","NC_3")
DEGs <- read.csv("2-runDEG_DESeq2_results.csv",header = T)
DEGs <- DEGs %>% filter(!group %in% "Stable")
DEG_up <- DEGs %>% filter(group %in% "Up") %>% top_n(100,log2FoldChange)
DEG_down <- DEGs %>% filter(group %in% "Down") %>% top_n(-100,log2FoldChange)
DEG_down %>% arrange(desc(abs(log2FoldChange))) %>% pull(X) %>% paste(.,collapse = ",")

DEG <- bind_rows(DEG_down,DEG_up)

group <- c(rep("Si",3),rep("Control",3))
names(group) <- colnames(TPM)
HeatmapAnnotation <- ComplexHeatmap::HeatmapAnnotation(
    Group = group,
    col = list(
        Group = c(
            "Si" = "#f31616",
            "Control" = "#0b82f1"  
        )
    ),
    gp = grid::gpar(col = "white"),
    show_annotation_name = FALSE,
    annotation_name_align = TRUE,
    annotation_name_side = "left"
)
ck_all <- clusterData(
    exp = TPM[DEG$X,],
    cluster.method = "mfuzz",
    cluster.num = 2,
    scaleData = TRUE
)

temp <- DEG
colnames(temp)[1] <- "gene" 
ck_all$wide.res <- left_join(ck_all$wide.res,temp %>% dplyr::select(gene,group),by = "gene")  %>% 
  mutate(cluster = case_when(group %in% "Down" ~ 1,
                             group %in% "Up" ~ 2)) %>% 
  dplyr::select(-group) %>% 
  arrange(cluster) 

ck_all$long.res <- ck_all$wide.res  %>% 
  pivot_longer(cols = c("NC_1","NC_2","NC_3","Si_1","Si_2","Si_3"),names_to = "cell_type",values_to = "norm_value") %>% 
  mutate(cluster_name = case_when(cluster == 1 ~ "cluster 1 (100)",
                                  cluster == 2 ~ "cluster 2 (100)")) 


jpeg("top100DEGs.jpg", height = 15, width = 20, units = "cm", res = 300)
visCluster(
    object = ck_all,
    plot.type = "heatmap",
    genes.gp = c("italic", 12, col = "orange"),
    markGenes.side = "right",
    show_row_names = FALSE,
    show_row_dend = FALSE,
    HeatmapAnnotation = HeatmapAnnotation,
    column.split = group,
    annnoblock.text = TRUE,
    add.mline = FALSE,
    panel.arg = c(2, 0.25, 4, "grey90", NA),
    lgd.label = names(table(group)))
dev.off()


enrichKEGG <- enrichCluster(
    object = ck_all,
    OrgDb = org.Hs.eg.db,
    type = "KEGG",
    pvalueCutoff = 0.05,
    topn = 100,
    seed = 111111
)
enrichBP <- enrichCluster(
    object = ck_all,
    OrgDb = org.Hs.eg.db,
    type = "BP",
    pvalueCutoff = 0.05,
    topn = 100,
    seed = 111111
)
enrichMF <- enrichCluster(
    object = ck_all,
    OrgDb = org.Hs.eg.db,
    type = "MF",
    pvalueCutoff = 0.05,
    topn = 100,
    seed = 111111
)
enrichCC <- enrichCluster(
    object = ck_all,
    OrgDb = org.Hs.eg.db,
    type = "CC",
    pvalueCutoff = 0.05,
    topn = 100,
    seed = 111111
)

markGenes <- c()
pdf("top100DEGs_enrichment.pdf", height =5, width = 12)
visCluster(
    object = ck_all,
    plot.type = "both",
    markGenes = markGenes,
    genes.gp = c("italic", 10, col = "orange"),
    markGenes.side = "right",
    show_row_names = FALSE,
    show_row_dend = FALSE,
    HeatmapAnnotation = HeatmapAnnotation,
    sample.order = c("Si_1","Si_2","Si_3","NC_1","NC_2","NC_3"),
    annnoblock.text = TRUE,
    panel.arg = c(2, 0.25, 4, "grey90", NA),
    annoTerm.data = enrichCC %>% group_by(group) %>% top_n(10, dplyr::desc(pvalue)),
    go.col = rep(c("#ff0000","#1F77B4FF","#ff0000","#1F77B4FF","#FF7F0EFF"), c(2,1,1,6,3)),
    lgd.label = names(table(group)),
    ctAnno.col = ggsci::pal_d3()(3),
    term.text.limit = c(10, 18),
    annoKegg.data = enrichKEGG %>% group_by(group) %>% top_n(10, dplyr::desc(pvalue)),
    annoTerm.mside = "right",
    kegg.col = rep(c("#1F77B4FF","#ff0000","#1F77B4FF","#FF7F0EFF"), c(5,3,2,10)),
    by.kegg = "anno_link",
    by.go = "anno_link",
    add.line = TRUE,
    line.side = "none"
)
dev.off()

