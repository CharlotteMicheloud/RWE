###################################################################################
#				                     	  IBKF Muenster                                     #
###################################################################################
#                                                                                 #
# Projekt/Study: The Sceptical $P$-value -- A Review                              #
# Ersterstellung (Author, Date, Institute): Koeppe, J., 06.05.24, IBKF            #
# Purpose of the programm: creating forest plot for replicated studiesa           #
###################################################################################

###################################################################################
# packages
#devtools::install_github("wilkelab/ungeviz")

library(tidyverse)
#library(openxlsx)
library(moments)
library(ggplot2)
library(ungeviz)
library(patchwork)
library(ggpubr)

#install.packages("ungeviz")


###################################################################################
##1. data preparation (code from  R. Heyard)

## Reading data through two separate CSVs: one with the results (one row per
## result, several rows per trial) and one CSV with the covariates (one row per
## trial)

#data path
path = "" #has to be defined.

setwd(path)

#load data

data <- read.csv2("all_trials_results.csv", header = TRUE) %>%
  group_by(Trial, Type) %>%
  mutate(Estimate = as.numeric(str_trim(str_split(pe_95ci,
                                                  "\\(", simplify = TRUE)[1])),
         Lower = as.numeric(
           str_split(str_split(pe_95ci, "\\(", simplify = TRUE)[2],
                     ",", simplify = TRUE)[1]),
         Upper = as.numeric(
           str_remove(str_trim(str_split(str_split(pe_95ci,
                                                   "\\(", simplify = TRUE)[2],
                                         ",", simplify = TRUE)[2]), "\\)"))) %>%
  ungroup() %>%
  mutate(Type =
           case_when(Type == "Pooled primary database study result" ~ "Pooled RWE",
                     Type == "Pooled intention-to-treat analysis" ~ "PooledITT",
                     Type == "RCT result" ~ "RCT",
                     TRUE ~ Type)) %>%
  select(-pe_95ci) %>%
  left_join(read.csv2("all_trials_covariates.csv", header = TRUE),
            by = c("Trial")) %>%
  rename(study_type = type) %>%
  select(Trial,Type,Estimate, Lower, Upper, study_type)

## Get the margins to compare the results to:
margins <- read.csv2("margins.csv", header = TRUE)
data <- data %>%
  left_join(margins,
            by = "Trial") %>%
  mutate(margin = as.numeric(margin))

# Delete those rows with missing Estimates:
data <- data %>%
  filter(!is.na(Estimate))

# ISAR-REACT5 and VERO do not have a result for RWE pooled, because a chi-square
# test indicated that results were heterogeneous across databases.
 data <- data %>%
   filter(!Trial %in% c("ISAR-REACT5", "VERO"))


#only data for RCT and Pooled. Further results are deleted
 data <- data %>%
   filter(Type != "PooledITT")

 data <- data %>%
   filter(Type != "Optum")

 data <- data %>%
   filter(Type != "MarketScan")

 data <- data %>%
   filter(Type != "Medicare")

#######################################################################################################
######separate plots for SUP and NI study types
####Sup:

data_Sup<-data[data$study_type == "Sup",]

data_Sup <- data_Sup[order(data_Sup$Trial),]

#adding line-number to data (for Forest-plot)
line_num<-seq(1,nrow(data_Sup))
data_Sup <- data_Sup %>% cbind(line_num)
data_Sup[data_Sup$Type== "Pooled RWE",]$line_num= data_Sup[data_Sup$Type== "Pooled RWE",]$line_num-0.5

studyNames <- unique(data_Sup$Trial) #for labels
as.vector(data_Sup[data_Sup$Type == "RCT",]$line_num)   #for labels

 p_Sup =
   ggplot(data=data_Sup, aes(x =(-line_num),y = Estimate, ymin = Lower, ymax =Upper) )+
   geom_hline(aes(),yintercept =1, linetype=2)+
   geom_pointrange(aes(col=Type))+
   coord_flip()+
   geom_linerange(ymin =-Inf, ymax =Inf, color = 'gray', size=0.3) +
   scale_color_manual(values = c("darkblue","black"),name="Type of study",na.translate = F)+
   xlab(' ')+ ylab("HR with 95%CI (logarithmic scale)")+
   geom_errorbar(aes(ymin= Lower, ymax= Upper,col=Type),width=0.3,cex=0.6)+
   scale_x_continuous(breaks= c(-1,-3,-5,-7,-9,-11,-13,-15,-17,-19,-21),
         labels= c(
          "-1"  = "DAPA-CKD",
          "-3"  = "HORIZON-PIVOTAL",
          "-5"  = "IMPACT",
          "-7"  = "INSPIRE",
          "-9"  = "P04334",
          "-11" = "PARADIGM-HF",
          "-13" = "PLATO",
          "-15" = "POET-COPD",
          "-17" = "PRONOUNCE",
          "-19" = "TRANSCEND",
          "-21" = "TRITON-TIMI"
         ))+
   scale_y_continuous(breaks=c(0.1, 0.25, 0.5, 1,2,4), trans = "log10")+
   theme(panel.grid.major.y = element_blank(),
         panel.grid.major.x = element_blank()) +
   theme(axis.text.y = element_text(face="bold", size=10),
         axis.ticks.y = element_line(size=0, colour = 'white'),
         panel.background = element_blank(),
         panel.border = element_rect(colour='black',fill=NA, size = 0.5),
         legend.background = element_rect(linetype = 1, size = 0.3, colour = 1),
   )+
   theme(legend.key = element_rect(fill = "white", color = NA),
         legend.position = c(0.8,0.9))

#Non-Inf separated with margin
data_NI<-data[data$study_type == "NI",]
data_NI<-data_NI[order(data_NI$Trial),]

#adding line-number to data (for Forest-plot)
line_num<-seq(1,nrow(data_NI))
data_NI <- data_NI %>% cbind(line_num)
data_NI[data_NI$Type== "Pooled RWE",]$line_num= data_NI[data_NI$Type== "Pooled RWE",]$line_num-0.5

studyNames <- unique(data_NI$Trial) #for labels
as.vector(data_NI[data_NI$Type == "RCT",]$line_num)   #for labels

p_NI =
  ggplot(data=data_NI, aes(x =(-line_num),y = Estimate, ymin = Lower, ymax =Upper) )+
  geom_pointrange(aes(col=Type))+
  coord_flip()+
  geom_linerange(ymin =-Inf, ymax =Inf, color = 'gray', size=0.3) +
  #adding margins
  geom_segment(
    data=tibble(),
    aes(
      x=-data_NI[data_NI$Type== "Pooled RWE",]$line_num+1.,
      y=data_NI[data_NI$Type== "Pooled RWE",]$margin,
      xend=-data_NI[data_NI$Type== "Pooled RWE",]$line_num-0.5,
      yend=data_NI[data_NI$Type== "Pooled RWE",]$margin
    ),
    inherit.aes=FALSE,linetype=2
  )+
  scale_color_manual(values = c("darkblue","black"),name=" ",na.translate = F)+
  xlab(' ')+ ylab("HR with 95%CI (logarithmic scale)")+
  geom_errorbar(aes(ymin= Lower, ymax= Upper,col=Type),width=0.3,cex=0.6)+
  scale_x_continuous(breaks= c(-1,-3,-5,-7,-9,-11,-13,-15,-17,-19,-21,-23,-25,-27,-29,-31,-33,-35),
                     labels= c(
                       "-1"  = "AMPLIFY",
                       "-3"  = "ARISTOTLE",
                       "-5"  = "CANVAS",
                       "-7"  = "CARMELINA",
                       "-9"  = "CAROLINA",
                       "-11" = "D5896",
                       "-13" = "DECLARE",
                       "-15" = "EINSTEIN-DVT",
                       "-17" = "EINSTEIN-PE",
                       "-19" = "EMPA-REG",
                       "-21" = "LEADER",
                       "-23" = "ON-TARGET",
                       "-25" = "RE-LY",
                       "-27" = "RECORD1",
                       "-29" = "RECOVER II",
                       "-31" = "ROCKET-AF",
                       "-33" = "SAVOR-TIMI",
                       "-35" = "TECOS"
                     ))+
  scale_y_continuous(breaks=c(0.1, 0.25, 0.5, 1,2,4), trans = "log10")+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank()) +
  theme(axis.text.y = element_text(face="bold", size=10),
        axis.ticks.y = element_line(size=0, colour = 'white'),
        panel.background = element_blank(),
        panel.border = element_rect(colour='black',fill=NA, size = 0.5),
        legend.background = element_rect(linetype = 1, size = 0.3, colour = 1),
  )+
  theme(legend.key = element_rect(fill = "white", color = NA),
        legend.position = c(0.12,0.9))


#both plots in one pdf:

layout <- "
AABB
AABB
AABB
"
p_total<-p_Sup +  labs(title = "Superiority trials") +
  p_NI +
  theme(legend.position="none")+  labs(title = "Noninferiority trials") +
  plot_layout(design=layout)


#output path
path_out = "" #has to be defined.


#Save plot:
ggsave('2024_05_06_Forest_pooled_both.pdf', plot = p_total, device = NULL,
       path = path_out,
       scale = 1,width =25, height = 20, units = "cm")


