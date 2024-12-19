
library(tidyverse)
library(ggplot2)


# principal component analysis --------------------------------------------


eigenvec <- read.table("pca.eigevec", quote="\"", comment.char="")
colnames(eigenvec)<-c("FID","sample",paste0("PC",1:10))
#eigenvec <- eigenvec[-c(25,38),]
write.table(eigenvec[2:ncol(eigenvec)],file = "pca.eigenvector.txt",sep = "\t",row.names = F,col.names = T,quote = F)

eigenval <- read.table("pca.eigenval", quote="\"", comment.char="")
pcs<-paste0("pc",1:nrow(eigenval))

percentage<-eigenval$V1/sum(eigenval$V1)*100
eigenval_df<-as.data.frame(cbind(pcs,eigenval[,1],percentage),stringsAsFactors = F)
names(eigenval_df)<-c("pcs","variance","proportation")
eigenval_df$variance<-as.numeric(eigenval_df$variance)
eigenval_df$proportation<-as.numeric(eigenval_df$proportation)

group_info <- read.csv("clinical_phenotype")

group_info <- group_info %>% mutate(color=if_else(outcome==1,"#9ACD32","#FF4500"),
                      pch=if_else(outcome==1,15,16),outcome=if_else(outcome==1,"control","case"))
#group_info$outcome <- as.factor(group_info$outcome)
#group_info <- group_info[-c(25,27),]
pdf("pca_20230322.pdf", width = 8,height = 6)
ggplot(data = eigenvec, aes(x = PC1, y = PC2, group = group_info$outcome)) +
  geom_point(alpha = 1,col=group_info$color,pch=group_info$pch)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1.5))+ #边框
  theme(legend.text = element_text(size = 16),legend.title = element_blank())+ #图例
  xlab(paste0("PC1 (", round(eigenval_df[eigenval_df$pcs == "pc1", 3], 2), "%)"))+ 
  ylab(paste0("PC2 (", round(eigenval_df[eigenval_df$pcs == "pc2", 3], 2), "%)")) +
  theme(axis.title.x = element_text(face = "bold", size = 18, colour = "black"),
        axis.title.y = element_text(face = "bold", size = 18, colour = "black"),
        axis.text.x = element_text(size = 14,face = "bold", colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black"))+
  stat_ellipse(aes(fill = group_info$outcome), geom = 'polygon', level = 0.99, alpha = 0.1, show.legend = T) +
  scale_fill_manual(values = c('orange', 'purple'))
dev.off()


##significance test of PCA
data <- merge(eigenvec,group_info,by = "sample")

data <- data %>% mutate(outcome=if_else(outcome=="case",2,1))
summary(glm(data$outcome ~ data$PC5))
