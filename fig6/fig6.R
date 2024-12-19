library(tableone)
library(xlsx)
library(tidyverse)
library(mice)


# datainput_and_cleaning --------------------------------------------------

raw_data <- read.xlsx("data.fig6.xlsx",header = T,sheetIndex="with_genotyping") %>% as_tibble() %>% 
  select(-诊断) %>%
  rename("death"="结局.存活1.死亡0.","FT3"="游离三碘甲腺原氨酸.FT3.","FT4"="游离甲状腺素.FT4.",
         "PTH"="甲状旁腺激素.PTH.","T3"="三碘甲腺原氨酸.T3.","T4"="甲状腺素.T4.","TG"="甲状腺球蛋白.TG.")

factor_var <- c("carrier","SIB","性别","病因","严重度","death","ARDS","AKI","胰腺坏死组织感染","腹腔高压","肝功能不全")
quant_var <- setdiff(colnames(raw_data)[-c(1,2)],factor_var)

clean_data <- raw_data %>% 
  mutate_all(na_if,"")  %>% 
  mutate(carrier=case_when(carrier==0~"not_carrier",
                           carrier==1~"carrier")) %>% 
  mutate(性别=case_when(性别=="男性"~"male",
                          性别=="女性"~"female")) %>%
  mutate(病因=case_when(病因==1~"HTG",
                        病因==2~"Alcohol",
                        病因==3~"biliary",
                        病因==4~"others")) %>%
  mutate(严重度=case_when(严重度==1~"mild",
                          严重度==2~"moderate",
                          严重度==3~"severe",
                          严重度==4~"critical")) %>%
  mutate(death=case_when(death==1~"survive",
                           death==0~"dead")) %>% 
  mutate(across(.cols = all_of(factor_var),as_factor)) %>% 
  mutate(across(.cols = all_of(quant_var),as.numeric))
  

glimpse(clean_data)


# comparison_and_create_table1 --------------------------------------------


tab1 <- CreateTableOne(vars = colnames(clean_data)[-c(1,2,3)], factorVars = factor_var[-1], 
                       data = clean_data, strata = "carrier", addOverall = TRUE)
tab1 <- print(tab1,showAllLevels=T,quote = F, noSpaces = TRUE,nonnormal=colnames(clean_data)[-c(1,2,3)] )

tab1 <- tab1  %>% as.data.frame()  %>% rownames_to_column(var = "rowname")
write_excel_csv(tab1,file = "./manuscripts/carrier_vs_noncarrier.csv")


# correlation_analysis_and_visualization -----------------------------------


library(psych)

sig_var <- c("SIB","性别","严重度","death","AKI","胰腺坏死组织感染","凝血酶时间","平均血红蛋白浓度","红细胞计数",
             "红细胞分布宽度","血红蛋白","红细胞比容","白细胞计数","谷氨酰转肽酶","尿素",
             "胱抑素C","钠","钙")

corr_analysis_data <- cbind(clean_data[,"carrier"],clean_data[,sig_var])

png('correlation_of_variables_carrier.png',units = "cm",res=300,width =100,height = 100)
pairs.panels(corr_analysis_data,scale = T,smooth = T,method = "kendall",stars = T)
dev.off()




# visualization_with_box-plot -----------------------------------------------

library(hrbrthemes)
library(viridis)
library(ggpubr)

corr_analysis_data <- cbind(clean_data[,"carrier"],clean_data[,sig_var])
p1 <-   corr_analysis_data %>% rename(value=sig_var[18]) %>% 
  #   filter(value<500) %>% 
  ggplot(aes(x=as_factor(carrier), y=value, fill=carrier)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  ggtitle(sig_var[18]) +
  xlab("")+ylab("")+
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  stat_compare_means(method="t.test",label = "p.format",label.x.npc = "center")+
  stat_n_text()+
  guides(fill="none")


png(paste0("carrier_info/",sig_var[18],".png"),units = "cm",res=600,width =20,height =15)
p1
dev.off() 

png(paste0("carrier_info/","SIB",".png"),units = "cm",res=600,width =20,height =15)
corr_analysis_data %>% mutate(SIB=case_when(SIB==0~"No",SIB==1~"Yes")) %>% 
  ggplot( aes(x = SIB,fill=SIB)) +
  geom_bar( )+
  theme_bw()+
  facet_grid(.~carrier)+
  guides(fill="none")
dev.off()   

png(paste0("carrier_info/","性别",".png"),units = "cm",res=600,width =20,height =15)
corr_analysis_data %>% 
  ggplot( aes(x = 性别,fill=性别)) +
  geom_bar( )+
  theme_bw()+
  facet_grid(.~carrier)+
  guides(fill="none")
dev.off()   

png(paste0("carrier_info/","严重度",".png"),units = "cm",res=600,width =20,height =15)
corr_analysis_data %>% 
  ggplot( aes(x = 严重度,fill=严重度)) +
  geom_bar( )+
  theme_bw()+
  facet_grid(.~carrier)+
  guides(fill="none")
dev.off()  

png(paste0("carrier_info/","death",".png"),units = "cm",res=600,width =20,height =15)
corr_analysis_data %>% 
  ggplot( aes(x = death,fill=death)) +
  geom_bar( )+
  theme_bw()+
  facet_grid(.~carrier)+
  guides(fill="none")
dev.off()   

png(paste0("carrier_info/","AKI",".png"),units = "cm",res=600,width =20,height =15)
corr_analysis_data %>% filter(!is.na(AKI)) %>% 
  ggplot( aes(x = AKI,fill=AKI)) +
  geom_bar( )+
  theme_bw()+
  facet_grid(.~carrier)+
  guides(fill="none")
dev.off()   

png(paste0("carrier_info/","胰腺坏死组织感染",".png"),units = "cm",res=600,width =20,height =15)
corr_analysis_data %>% filter(!is.na(胰腺坏死组织感染)) %>%
  ggplot( aes(x = 胰腺坏死组织感染,fill=胰腺坏死组织感染)) +
  geom_bar( )+
  theme_bw()+
  facet_grid(.~carrier)+
  guides(fill="none")
dev.off()   




# visualization_costumed_clinical_deaths---------------------------------------------------
 
#remotes::install_github("kcuilla/reactablefmtr")
library("reactablefmtr")

sig_var_dichotomous <- c("Spontaneous intraabdominal bleeding","Critical severity","Death",
                         "Acute kidney injury","Infected pancreatic necrosis")

carrier_pct <- c(0.25,0.431,0.118,0.465,0.562)
non_carrier_pct <- c(0.044,0.202,0.044,0.335,0.269)
p_value <- c("<0.001","<0.001","0.021","0.023","<0.001") %>% as.character()

summary <- data.frame(deaths=sig_var_dichotomous,carrier=carrier_pct,non_carrier=non_carrier_pct,p_value)

cli_death <- reactable(summary,
          pagination = F,
          highlight = T,
          theme = fivethirtyeight(),
          defaultSortOrder = "desc",
          style = list(fontFamily = "Arial"),
          defaultColDef = colDef(maxWidth = 150,format=colFormat(digits = 0)),
          columns = list(deaths=colDef(maxWidth = 300),
                         carrier = colDef(cell = data_bars(summary,  fill_color = "#E64B35cc",
                                                               align_bars = "right", 
                                                               text_position = "inside-end", 
                                                               box_shadow = F,round_edges = F,
                                                               number_fmt = scales::percent,
                                                           background = "transparent"
                                                               ),
                                          maxWidth = 200,
                                          align = "center"
                                              ),
                         non_carrier = colDef(cell = data_bars(summary,  fill_color = "#4DBBD5CC",
                                                                   align_bars = "left",
                                                                   text_position = "inside-end",
                                                                   box_shadow = F,round_edges = F,
                                                                   number_fmt = scales::percent,
                                                               background = "transparent"
                                                                   ),
                                              maxWidth = 200,
                                              align = "center"
                                                  ),
                       p_value=colDef(maxWidth = 70,align = "center")
                         )
          )
cli_death
save_reactable(cli_death, "carrier_info/clin_death.png")
save_reactable(cli_death, "carrier_info/clin_death.html")


# visualization_costumed_lab_tests_part1 ----------------------------------


sig_var_quantities <- c("Thrombin time (second)","Red blood cell count (×10^12/L)",
                         "Red blood Cell distribution width (%)","Hematocrit (%)","White blood cell count (×10^9/L)",
                        "Blood urea nitrogen (mmol/L)","Cystatin C (mg/L)","Phosphorus (mmol/L)","Alcium (mmol/L)")

carrier_pct <- c(19.24,3.22,15.07,29,10.47,8.95,1.30,1.18,2.05)
non_carrier_pct <- c(16.97,3.57,14.23,33,12.31,6.44,0.87,0.98,1.98)
p_value <- c("0.014","<0.001","<0.001","0.001","0.007","0.004","<0.001","0.001","0.031") %>% as.character()

summary <- data.frame(lab_tests=sig_var_quantities,carrier=carrier_pct,non_carrier=non_carrier_pct,p_value)

lab_test <- reactable(summary,
                         pagination = FALSE,
                         highlight = T,
                         theme = fivethirtyeight(),
                         style = list(fontFamily = "Arial"),
                         defaultSortOrder = "desc",
                         defaultColDef = colDef(maxWidth = 150,format=colFormat(digits = 0)),
                         columns = list(lab_tests=colDef(maxWidth = 300),
                                        carrier = colDef(cell = data_bars(summary,  fill_color = "#E64B35cc",
                                                                          align_bars = "right", 
                                                                          text_position = "inside-end", 
                                                                          force_outside = c(0,2.5),
                                                                          box_shadow = F,round_edges = F,
                                                                          max_value = 15,
                                                                          background = "transparent"
                                        ),
                                        maxWidth = 200,
                                        align = "center"
                                        ),
                                        non_carrier = colDef(cell = data_bars(summary,  fill_color = "#4DBBD5CC",
                                                                              align_bars = "left",
                                                                              text_position = "inside-end",
                                                                              force_outside = c(0,2.5),
                                                                              box_shadow = F,round_edges = F,
                                                                              max_value = 15,
                                                                              background = "transparent"
                                        ),
                                        maxWidth = 200,
                                        align = "center"
                                        ),
                                        p_value=colDef(maxWidth = 70,align = "center")
                         )
)
lab_test
save_reactable(lab_test, "./manuscripts/carrier_info/clin_death_part1.png")


# visualization_costumed_lab_tests_part2 ----------------------------------



sig_var_quantities <- c("Mean corpuscular hemoglobin concentration (g/L)","Hemoglobin (g/L)","Glutamyl transpeptidase (U/L)")

carrier_pct <- c(325.72,95.36,125.95)
non_carrier_pct <- c(332.38,107.30,84.64)
p_value <- c("0.007","<0.001","0.013") %>% as.character()

summary <- data.frame(lab_tests=sig_var_quantities,carrier=carrier_pct,non_carrier=non_carrier_pct,p_value)

lab_test <- reactable(summary,
                      pagination = FALSE,
                      highlight = T,
                      theme=fivethirtyeight(),
                      defaultSortOrder = "desc",
                      style = list(fontFamily = "Arial"),
                      defaultColDef = colDef(maxWidth = 150,format=colFormat(digits = 0)),
                      columns = list(lab_tests=colDef(maxWidth = 300),
                                     carrier = colDef(cell = data_bars(summary,  fill_color = "#E64B35cc",
                                                                       align_bars = "right", 
                                                                       text_position = "inside-end", 
                                                                       force_outside = c(0,2.5),
                                                                       box_shadow = F,round_edges = F,
                                                                       max_value = 130,
                                                                       background = "transparent"
                                     ),
                                     maxWidth = 200,
                                     align = "center"
                                     ),
                                     non_carrier = colDef(cell = data_bars(summary,  fill_color = "#4DBBD5CC",
                                                                           align_bars = "left",
                                                                           text_position = "inside-end",
                                                                           force_outside = c(0,2.5),
                                                                           box_shadow = F,round_edges = F,
                                                                           max_value = 130,
                                                                           background = "transparent"
                                     ),
                                     maxWidth = 200,
                                     align = "center"
                                     ),
                                     p_value=colDef(maxWidth = 70,align = "center")
                      )
)
lab_test
save_reactable(lab_test, "carrier_info/clin_death_part2.png")



# logistic_regression_and_visualization --------------------------

clean_data_backup <- clean_data



clean_data <- clean_data_backup
clean_data <- clean_data %>% 
  mutate(严重度=case_when(严重度=="mild"~0,
                          严重度=="moderate"~0,
                          严重度=="severe"~0,
                          严重度=="critical"~1)) %>%
  mutate(death=case_when(death=="survive"~0,
                         death=="dead"~1)) %>%
  mutate(across(.cols = 严重度,as_factor)) %>% 
  mutate(across(.cols = death,as_factor))


#SIB  严重度  death  ARDS  AKI  胰腺坏死组织感染  腹腔高压  肝功能不全
consequence <- "SIB"
forestplot_title <- "Spontaneous intraabdominal hemorrhage"

form1 <- as.formula(paste(consequence,"~","carrier"))
form2 <- as.formula(paste(consequence,"~","carrier+age"))
form3 <- as.formula(paste(consequence,"~","carrier+age+性别"))
form4 <- as.formula(paste(consequence,"~","carrier+age+性别+病因"))
form5 <- as.formula(paste(consequence,"~","carrier+age+性别+病因+BMI"))
form6 <- as.formula(paste(consequence,"~","age+性别+病因"))

mod1 <- glm(form1,data = clean_data,family = "binomial",control = list(maxit = 100));summary(mod1)
mod2 <- glm(form2,data = clean_data,family = "binomial",control = list(maxit = 100));summary(mod2)
mod3 <- glm(form3,data = clean_data,family = "binomial",control = list(maxit = 100));summary(mod3)
mod4 <- glm(form4,data = clean_data,family = "binomial",control = list(maxit = 100));summary(mod4)
mod5 <- glm(form5,data = clean_data,family = "binomial",control = list(maxit = 100));summary(mod5)
mod6 <- glm(form6,data = clean_data,family = "binomial",control = list(maxit = 100));summary(mod6)


###forest_Plot
library(forestplot)

windowsFonts(Arial=windowsFont("Arial"))

data_for_plot <- 
  structure(list(
    mean  = c(exp(coef(mod1))["carriercarrier"],exp(coef(mod2))["carriercarrier"],exp(coef(mod3))["carriercarrier"],exp(coef(mod4))["carriercarrier"],exp(coef(mod5))["carriercarrier"]), 
    lower = c(exp(confint(mod1))[2,1],exp(confint(mod2))[2,1],exp(confint(mod3))[2,1],exp(confint(mod4))[2,1],exp(confint(mod5))[2,1]),
    upper = c(exp(confint(mod1))[2,2],exp(confint(mod2))[2,2],exp(confint(mod3))[2,2],exp(confint(mod4))[2,2],exp(confint(mod5))[2,2])),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -5L), 
    class = "data.frame")
data_for_plot

row_names <- list(list("No adjustment", "Adjusted with age","Adjusted with age/gender","Adjusted with age/gender/etiology","Adjusted with age/gender/etiology/BMI"))
f1 <- forestplot(labeltext=row_names,
                 graph.pos = 2,
                 txt_gp = fpTxtGp(label = gpar(fontfamily = "Arial",cex=2.5),
                                  ticks = gpar(fontfamily = "Arial", cex =2)),
                 hrzl_lines = gpar(col = "#444444"),
                 new_page = TRUE,
                 data_for_plot,
                 clip=c(0.1,40),
                 is.summary = c(rep(FALSE,4)),
                 xlog = TRUE, 
                 col = fpColors(box = "#BC3C28",
                                line = "black",
                                summary = "#BC3C28"),
                 vertices = T,
                 title = forestplot_title)
f1

png(paste0("./manuscripts/carrier_info/",forestplot_title,"_forest.png"),units = "cm",res=600,width = 27,height = 18)
f1
dev.off()

##ROC plot
pre4 <- predict(mod4,type='response')
pre6 <- predict(mod6,type='response')

library(pROC)

clean_data <- clean_data[-34,]

rocobj4<-roc(clean_data[[consequence]],pre4,smooth=F, ci=TRUE)
auc4<-round(auc(clean_data[[consequence]],pre4),4)
rocobj6<-roc(clean_data[[consequence]],pre6,smooth=F, ci=TRUE)
auc6<-round(auc(clean_data[[consequence]],pre6),4)

p <- roc.test(rocobj4,rocobj6,method = 'delong');p

rocobj4<-roc(clean_data[[consequence]],pre4,smooth=T, ci=TRUE)
auc4<-round(auc(clean_data[[consequence]],pre4),4)
rocobj6<-roc(clean_data[[consequence]],pre6,smooth=T, ci=TRUE)
auc6<-round(auc(clean_data[[consequence]],pre6),4)

list <- list(with_mutation_info=rocobj4,without_mutation_info=rocobj6) 

library(ggsci)
p1<- ggroc(list,legacy.axes = TRUE,size=1,alpha=1)+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
               color="darkgrey", linetype=4)+
  theme_bw()+
  ggtitle(forestplot_title)+
  ggsci::scale_color_lancet()+
  annotate("text",x=0.73,y=0.25,label=paste("with mutation info-AUC = ", round(auc4,digits = 3)))+
  annotate("text",x=0.74,y=0.18,label=paste("without mutation info-AUC = ", round(auc6,digits = 3)))+
  annotate("text",x=0.715,y=0.11,label=paste("delong'test p_value = ", round(p$p.value,digits = 3)))+
  theme(axis.title.x = element_text(size = 20, family = "myFont", 
                                    color = "black", face = "bold"),
        axis.title.y = element_text(size = 20, family = "myFont", 
                                    color = "black", face = "bold"),
        axis.text.x = element_text(size = 15, family = "myFont",
                                   color = "black", face = "bold"),
        axis.text.y = element_text(size = 15, family = "myFont",
                                   color = "black", face = "bold"),
        legend.position = c(0.8,0.5),
        legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5,size = 20,face = "bold"))
p1

png(paste0("carrier_info/",forestplot_title,"_roc.jpg"),units="in", width=8, height=5,res=600)
p1
dev.off()
