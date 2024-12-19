##frequency histogram

library(tidyverse)
library(hrbrthemes)
phred <- read.csv("data.3D.csv")
colnames(phred)[1] <- "phred"


p <- phred %>%
  ggplot( aes(x=phred)) +
  geom_histogram( binwidth=0.8, fill="#69b3a2", color="#e9ecef", alpha=0.9,) +
  ggtitle("") +
  theme_ipsum() +
  theme(
    plot.title = element_text(size=15)
  )
p


  

png("CADD_phred_density.png",units="in", width=8, height=5,res=600)
phred %>%
  ggplot( aes(x=phred)) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8)+
  theme_classic()+
  xlab("Phred-scale score") + ylab("Density") + ggtitle("")+
  theme(axis.title.x = element_text(size = 20, family = "myFont", 
                                    color = "black", face = "bold"),
        axis.title.y = element_text(size = 20, family = "myFont", 
                                    color = "black", face = "bold"),
        axis.text.x = element_text(size = 15, family = "myFont",
                                   color = "black", face = "bold"),
        axis.text.y = element_text(size = 15, family = "myFont",
                                   color = "black", face = "bold"))
dev.off()
