
library(tidyverse)
library(ggrepel)

raw_snps <- read.csv("data.3C.csv") %>% 
  unite("ID",Chr,Position,Ref,Alt,sep = ":",na.rm = FALSE, remove = FALSE) 

a <- replicate(22,"chr")  
b <- c(1:22) 
autosomal_chr <- paste0(a,b)  

raw_snps <- raw_snps %>% filter(Chr %in% autosomal_chr) %>% 
  separate(Chr,into = c("nonuse","CHR_number"),sep = "r") %>% 
  select(ID,CHR_number, Position,p_ad)%>% 
  rename(Chr=CHR_number) %>% 
  mutate(across(.cols=Chr,as.numeric)) %>% 
  arrange(Chr,Position) %>% 
  mutate( is_annotate=ifelse(-log10(p_ad)>8, "yes", "no"))


chr_len <- raw_snps %>% 
  group_by(Chr) %>% 
  summarise(chr_len=max(Position))

chr_pos <- chr_len  %>% 
  mutate(total = cumsum(as.numeric(chr_len)) - chr_len) %>%
  select(-chr_len) 

Snp_pos <- chr_pos %>%
  left_join(raw_snps, ., by="Chr") %>%
  arrange(Chr, Position) %>%
  mutate( BPcum = Position + total)


X_axis <-  Snp_pos %>% group_by(Chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

p <- ggplot(Snp_pos, aes(x=BPcum, y=-log10(p_ad))) +
  geom_point( aes(color=as.factor(Chr)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  scale_x_continuous( label = X_axis$Chr, breaks= X_axis$center ) +
  scale_y_continuous(expand = c(0, 0),limits = c(0,25) ) + 
  geom_hline(yintercept = c(6, -log10(0.05/nrow(Snp_pos))), color = c('green', 'red'), 
             linewidth = 1.2, linetype = c("dotted", "twodash")) +
  geom_point(data=subset(Snp_pos, is_annotate=="yes"), color="orange", size=2)+
  geom_label_repel( data=subset(Snp_pos, is_annotate=="yes"), aes(label=ID), size=2)+
  theme_bw() +
  theme(
    legend.position="none",
    panel.border = element_blank(),
    axis.line.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

p
png("manhattan.png",units="cm",res=1200,width =24,height =18)
p
dev.off() 
