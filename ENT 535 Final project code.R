library(tidyverse)
library(svglite)
library(Cairo)
library(ggplot2)
install.packages("rstatix")
library(rstatix)
library(dplyr)
library(ggpubr)

#I want to test if there is a significant difference in silica accumulation between WT and Bdgt43b2 

#This requires a two sample t test because I am testing differences between two independent groups

#### Import and clean up data file ####
#Import csv file
silica <- read.csv("Anderson Lab - Alisa Chernikova - Anderson, Charles T's files/2024 NSF GRFP/2024.10.11_BrachyLeaf_SilicaCellQuantification/2024.10.11_BrachySilicaCellQuantification_CorrectedTotalCellFluorescence.csv")

#Subset columns of interest (Genotype, LeafSide, CorrectedTotalCellFluorescence) into a new object
silica2 <- subset(silica, select = c(Genotype, LeafSide, CorrectedTotalCellFluorescence))

#Convert the object into a dataframe
silica3 <- as.data.frame(silica2)

#Calculate outliers
silica_outlier<-
  silica3 %>%
  group_by(Genotype, LeafSide) %>% 
  identify_outliers("CorrectedTotalCellFluorescence") 
#These are genuine biological outliers (fluorescence was much higher than fluorescence in the average silica cell) so I am comfortable in removing the outliers from my dataset

#Make a new dataframe with outliers removed and replaced with NAs
#I used the 1.5*IQR rule since it is the standard outlier test
silica_cleaned<- silica3 %>% 
  group_by(Genotype, LeafSide) %>% 
  mutate(CorrectedTotalCellFluorescence_cleaned = ifelse(CorrectedTotalCellFluorescence>quantile(CorrectedTotalCellFluorescence, 0.75, na.rm = TRUE)+1.5*IQR(CorrectedTotalCellFluorescence, na.rm = TRUE) | CorrectedTotalCellFluorescence<quantile(CorrectedTotalCellFluorescence, 0.25, na.rm = TRUE)-1.5*IQR(CorrectedTotalCellFluorescence, na.rm = TRUE), NA, CorrectedTotalCellFluorescence))


#### Assumption time ####
#Two sample ttest assumptions
#1: data are continuous 
#TRUE

#2: data are in independent and randomly sampled
#TRUE

#3: data are normal
hist(silica_cleaned$CorrectedTotalCellFluorescence_cleaned)
#Data does NOT have a normal distribution

#4: data has equal variance
#Null: data has equal variance (p>0.05)
#Alternate: das unqual variance (p<0.05)
bartlett.test(CorrectedTotalCellFluorescence_cleaned~LeafSide,data=silica_cleaned)
#Bartlett's K-squared = 2.8131, df = 1, p-value = 0.0935
#p>0.05 so variances are equal


#### Statistical test time ####
#Since data are not normally distributed and do not have equal variance I will perform the Mann-Whitney test
test<-silica_cleaned %>% 
  group_by(LeafSide) %>% 
  wilcox_test(CorrectedTotalCellFluorescence_cleaned~Genotype) %>% 
  adjust_pvalue() %>%
  add_significance("p.adj") %>% 
  add_xy_position(x = "LeafSide", group = "Genotype", dodge = 1)

#Print Mann-Whitney test results
print(test)


#### Graph time ####
ggplot(silica_cleaned, aes(x=LeafSide, y=CorrectedTotalCellFluorescence_cleaned, color=Genotype))+
  geom_jitter(position=position_jitterdodge(jitter.width=0.1, dodge.width = 0.8), cex=2.5)+
  stat_summary(aes(group=Genotype), 
               fun=mean, 
               geom="crossbar", color="black", 
               position=position_dodge(width=0.8), width=0.4)+
  labs(y="Corrected Total Cell Fluorescence (a.u.)", face = "bold")+
  theme_light()+
  coord_cartesian(clip = "off")+
  scale_color_manual(values=c("olivedrab3", "forestgreen"))+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color="black", linewidth=0.5),
        axis.title.y = element_text(face="bold", 
                                    size = 13),
        axis.ticks = element_line(color = "black",
                                  linewidth = 1.5),
        axis.ticks.length = unit(0.3, "cm"),
        axis.text.y = element_text(face = "bold", 
                                   size = 13, 
                                   color="black"),
        axis.text.x = element_text(face = "bold", 
                                   size=13, 
                                   color="black",
                                   angle = 45, 
                                   hjust = 1),
        legend.title = element_text(face = "bold", 
                                    size=13, 
                                    color="black"),
        legend.text = element_text(face = "bold", 
                                   size=13, 
                                   color="black"),
        axis.title.x=element_blank(),
        panel.border = element_rect(color= "black", fill=NA, linewidth=1.5))+
  stat_pvalue_manual(test, step.group.by = "LeafSide", label = "p.adj.signif", step.increase = 2, tip.length = 0.02, fontface = "bold", label.size = 5, bracket.size = 1.1, bracket.nudge.y = 0.1e+06)+
  scale_x_discrete(labels=c("AB"="Abaxial", "AD"="Adaxial"))+
  scale_y_continuous(limits=c(0,5e+06), expand = c(0, 0))

#save the plot as an object in order to export the plot from R
silica_plot <- ggplot(silica_cleaned, aes(x=LeafSide, y=CorrectedTotalCellFluorescence_cleaned, color=Genotype))+
  geom_jitter(position=position_jitterdodge(jitter.width=0.1, dodge.width = 0.8), cex=2.5)+
  stat_summary(aes(group=Genotype), 
               fun=mean, 
               geom="crossbar", color="black", 
               position=position_dodge(width=0.8), width=0.4)+
  labs(y="Corrected Total Cell Fluorescence (a.u.)", face = "bold")+
  theme_light()+
  coord_cartesian(clip = "off")+
  scale_color_manual(values=c("olivedrab3", "forestgreen"))+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color="black", linewidth=0.5),
        axis.title.y = element_text(face="bold", 
                                    size = 13),
        axis.ticks = element_line(color = "black",
                                  linewidth = 1.5),
        axis.ticks.length = unit(0.3, "cm"),
        axis.text.y = element_text(face = "bold", 
                                   size = 13, 
                                   color="black"),
        axis.text.x = element_text(face = "bold", 
                                   size=13, 
                                   color="black",
                                   angle = 45, 
                                   hjust = 1),
        legend.title = element_text(face = "bold", 
                                    size=13, 
                                    color="black"),
        legend.text = element_text(face = "bold", 
                                   size=13, 
                                   color="black"),
        axis.title.x=element_blank(),
        panel.border = element_rect(color= "black", fill=NA, linewidth=1.5))+
  stat_pvalue_manual(test, step.group.by = "LeafSide", label = "p.adj.signif", step.increase = 2, tip.length = 0.02, fontface = "bold", label.size = 5, bracket.size = 1.1, bracket.nudge.y = 0.1e+06)+
  scale_x_discrete(labels=c("AB"="Abaxial", "AD"="Adaxial"))+
  scale_y_continuous(limits=c(0,5e+06), expand = c(0, 0))

#to save the file as a .png or .svg change the file ending in paste0()
ggsave(paste0("silica_plot.png"), plot=silica_plot, path="~/Desktop", dpi = 1000, width = 6, height = 4.16, units = "in")
