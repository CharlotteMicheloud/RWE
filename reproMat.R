# Reproducibility Material for the paper
# 'Assessing the replicability of RCTs in RWE emulations'
# from J. Köppe, C. Micheloud, S. Erdmann, R. Heyard and L. Held


library(ReplicationSuccess)
library(biostatUZH)
library(tidyverse)
library(moments)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(gplots)
library(moments)
library(patchwork)


# The three datasets used come from R. Heyard's GitLab page: 
# https://gitlab.com/heyardr/hte-in-rwe

###################################################################################
##1. data preparation (code from  R. Heyard)

## Reading data through two separate CSVs: one with the results (one row per
## result, several rows per trial) and one CSV with the covariates (one row per
## trial)

#load data
# setwd(path)
data <- read.csv2("data/all_trials_results.csv", header = TRUE) %>%
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
  left_join(read.csv2("data/all_trials_covariates.csv", header = TRUE),
            by = c("Trial")) %>%
  rename(study_type = type) %>%
  select(Trial,Type,Estimate, Lower, Upper, study_type)

## Get the margins to compare the results to:
margins <- read.csv2("data/margins.csv", header = TRUE)
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

studyNames <- unique(data$Trial)


################################################################################
######################  Figure 1: Forest plot  #################################
################################################################################



#######################################################################################################

# Superiority
data_Sup<-data[data$study_type == "Sup",]
data_Sup <- data_Sup[order(data_Sup$Trial),] #alphabetical order

# adding line-number to data (for Forest-plot)
line_num<-seq(1,nrow(data_Sup))
data_Sup <- data_Sup %>% cbind(line_num)
data_Sup[data_Sup$Type== "Pooled RWE",]$line_num= data_Sup[data_Sup$Type== "Pooled RWE",]$line_num-0.5

studyNames2 <- unique(data_Sup$Trial)                    #for labels
as.vector(data_Sup[data_Sup$Type == "RCT",]$line_num)   #for labels -> added as labels for ticks!

  ggplot(data=data_Sup, aes(x =(-line_num),y = Estimate, ymin = Lower, ymax =Upper) )+
  geom_hline(aes(),yintercept =1, linetype=2)+
  geom_pointrange(aes(col=Type))+
  coord_flip()+
  geom_linerange(ymin =-Inf, ymax =Inf, color = 'gray', size=0.3) +
  scale_color_manual(values = c("darkblue","black"),name="Type of study",na.translate = F)+
  xlab(' ')+ ylab("HR with 95%CI (logarithmic scale)") +
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


# Non-Inferiority 
  
  
data_NI<-data[data$study_type == "NI",]
data_NI<-data_NI[order(data_NI$Trial),]

#adding line-number to data (for Forest-plot)
line_num<-seq(1,nrow(data_NI))
data_NI <- data_NI %>% cbind(line_num)
data_NI[data_NI$Type== "Pooled RWE",]$line_num= data_NI[data_NI$Type== "Pooled RWE",]$line_num-0.5

studyNames3 <- unique(data_NI$Trial) #for labels
as.vector(data_NI[data_NI$Type == "RCT",]$line_num)   #for labels


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
  xlab(' ')+ ylab("HR with 95%CI (logarithmic scale)") +
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
  p_NI + theme(legend.position="none")+  labs(title = "Noninferiority trials") +
  plot_layout(design=layout)


################################################################################
######################  Figure 2: Sceptical p  #################################
################################################################################

  # @ Leo, please add your code



################################################################################
######################  Figure 3: Cond vs Pred pow  ############################
################################################################################
  data <- data %>%
    mutate(se = ci2se(lower = Lower, upper = Upper, ratio = TRUE),
           z = (log(Estimate) - log(margin))/se)
  # summary(data)
  
  
  
  ## Extract the results from the RCTs and the (pooled) results from the RWE
  ## studies
  dataRCT <- data %>%
    filter(Type == "RCT")
  dataRWE <- data %>%
    filter(Type == "Pooled RWE")
  
  # ISAR-REACT5 and VERO do not have results for the pooled RWE, so we do not
  # include them here:
  # TODO: or delete them from the start
  # ids <- which(studyNames %in% c("ISAR-REACT5", "VERO"))
  
  ## Get the results matrix
  res <- lapply(studyNames, function(name) {
    # For each study, extract the RCT and the pooled RWE
    RCT <- dataRCT %>%
      filter(Trial == name)
    RWE <- dataRWE %>%
      filter(Trial == name)
    
    # Extract the elements of the results:
    #   * the effect sizes (theta_o and theta_r)
    #   * the golden and the controlled sceptical p-values
    #   * the p-value of the 2 trials rule (max of original and replication)
    #   * the variance ratio
    #   * the standard errors and the z values (of original and replication)
    r <- c(thetahato = log(RCT$Estimate),
           thetahatr = log(RWE$Estimate),
           st_d = (log(RCT$Estimate) -
                     log(RWE$Estimate))/sqrt(RCT$se**2 + RWE$se**2),
           # CM: to avoid problem with PRONOUNCE
           golden_pS = ifelse(RCT$z < 0, 
                              pSceptical(zo = RCT$z,
                                         zr = RWE$z,
                                         c = RCT$se^2/RWE$se^2,
                                         type = "golden"),
                              1 - pSceptical(zo = RCT$z,
                                             zr = RWE$z,
                                             c = RCT$se^2/RWE$se^2,
                                             type = "golden")),
           # CM: to avoid problem with PRONOUNCE
           controlled_pS = ifelse(RCT$z < 0, 
                                  pSceptical(zo = RCT$z,
                                             zr = RWE$z,
                                             c = RCT$se^2/RWE$se^2,
                                             type = "controlled"),
                                  1 - pSceptical(zo = RCT$z,
                                                 zr = RWE$z,
                                                 c = RCT$se^2/RWE$se^2,
                                                 type = "controlled")), 
           TTR_p = NA,
           varianceRatio = RCT$se^2/RWE$se^2,
           se_o = RCT$se, se_r = RWE$se,
           zo = RCT$z, zr = RWE$z,
           # For the Meta-analysis criterion
           theta_MA = sum(c(log(RCT$Estimate)/(RCT$se)^2,
                            log(RWE$Estimate)/(RWE$se)^2))/sum(c(1/RCT$se^2,
                                                                 1/RWE$se^2)),
           se_MA = 1/sqrt(sum(c(1/RCT$se^2, 1/RWE$se^2))),
           margin = RCT$margin)
    return(r)
  })
  
  # Bind list together into a dataframe
  res <- do.call("rbind", res) %>%
    as.data.frame() %>%
    # Compute the Meta-analysis criterion using the NI-margin
    mutate(z_MA = (log(margin) - theta_MA)/se_MA, #(theta_MA - log(margin))/se_MA,
           MA_p = z2p(z_MA, alternative = "one.sided")) %>%
    # Add the name of the study as first column
    mutate(Study = studyNames) %>%
    select(Study, everything(), -z_MA)
  
  res$po <- z2p(res$zo, alternative = "less")
  res$pr <- z2p(res$zr, alternative = "less")
  
  res$TTR_p <- pmax(res$po, res$pr)
  
  
alpha <-  0.025
  study <- subset(res, Study == "TRITON-TIMI")
  
  # CI for original study
  CIo <- study$thetahato + c(-1, 1)*qnorm(1-alpha)*study$se_o
  
  # CIr <- study$thetahatr + c(-1, 1)*qnorm(1-alpha)*study$se_r
  
  
  # cond power for original effect estimate
  cond_pow <- powerSignificance(zo = study$zo, c = study$se_o^2/study$se_r^2,
                                level = 0.025, designPrior = "conditional",
                                alternative = "one.sided")
  
  pred_pow <- powerSignificance(zo = study$zo, c = study$se_o^2/study$se_r^2,
                                level = 0.025, designPrior = "predictive",
                                alternative = "one.sided")
  
  # cond power for limits of CI
  cond_powCIl <- powerSignificance(zo = CIo[1]/study$se_o, c = study$se_o^2/study$se_r^2,
                                   level = 0.025, designPrior = "conditional",
                                   alternative = "one.sided")
  
  cond_powCIu <- powerSignificance(zo = CIo[2]/study$se_o, c = study$se_o^2/study$se_r^2,
                                   level = 0.025, designPrior = "conditional",
                                   alternative = "one.sided")
  
  delta_seq <- seq(-0.4, 0, by = 0.001)
  cond_pow_plot <- powerSignificance(zo = delta_seq/study$se_o,
                                     c = study$se_o^2/study$se_r^2,
                                     alternative = "one.sided",
                                     designPrior ="conditional")
  
  cond_pow_plot[delta_seq >= 0] <- NA
  
  
  
  
  # simulation for distribution of cond power
  thetao_sample <- rnorm(n = 100000, mean = study$thetahato, sd = study$se_o)
  condPower_sample <- powerSignificance(zo = thetao_sample/study$se_o, c = study$se_o^2/study$se_r^2,
                                        level = 0.025, designPrior = "conditional",
                                        alternative = "one.sided")
  
  
  
  par(las = 1)
  layout_mat <- matrix(c(2, 0,
                         1, 3),
                       nrow = 2,
                       byrow = T)
  
  layout(layout_mat, width = c(3, 3), height = c(2.5, 3),
         respect = TRUE)
  
  par(mar = c(4.1,4,0,0))
  
  plot(delta_seq, cond_pow_plot*100, type = "l", lwd = 2,
       xlab = "Original effect estimate",
       ylab = "Conditional power (%)",
       xlim = c(-0.4, 0),
       yaxt = "n",
       cex.axis = 0.7)
  
  segments(x0 = study$thetahato, x1 = study$thetahato, y0 = 0, y1 = cond_pow*100, col = "red",
           lty = 2)
  
  segments(x0 = -2, x1 = study$thetahato, y0 = cond_pow*100, y1 = cond_pow*100, col = "red",
           lty = 2)
  
  
  plotCI(x = study$thetahato,
         y = 0.5,
         li = CIo[1],
         ui = CIo[2],
         col=2, lwd=2, add=TRUE, err="x",
         barcol = "darkgreen")
  
  segments(x0 = CIo[1], x1 = CIo[1], y0 = 0, y1 = cond_powCIl*100, col = "darkgreen",
           lty = 2)
  
  segments(x0 = -2, x1 = CIo[1], y0 = cond_powCIl*100, y1 = cond_powCIl*100, col = "darkgreen",
           lty = 2)
  
  
  segments(x0 = CIo[2], x1 = CIo[2], y0 = 0, y1 = cond_powCIu*100, col = "darkgreen",
           lty = 2)
  
  segments(x0 = -2, x1 = CIo[2], y0 = cond_powCIu*100, y1 = cond_powCIu*100, col = "darkgreen",
           lty = 2)
  
  
  axis(2, at = cond_powCIl*100, labels = round(cond_powCIl*100, 1),
       col.axis = "darkgreen",
       col.ticks= "darkgreen",
       cex.axis = 0.6)
  axis(2, at = cond_powCIu*100, labels = round(cond_powCIu*100, 1),
       col.axis = "darkgreen",
       col.ticks= "darkgreen",
       cex.axis = 0.6)
  axis(2, at = cond_pow*100, labels = round(cond_pow*100, 1),
       col.axis = "red",
       col.ticks= "red",
       cex.axis = 0.6)
  axis(2, at = pred_pow*100, labels = round(pred_pow*100, 1),
       col.axis = "blue",
       col.ticks= "blue",
       cex.axis = 0.6)
  
  axis(2, at = c(0, 20, 40, 60, 80), labels = c(0, 20, 40, 60, 80),
       cex.axis = 0.6)
  
  
  plotCI(x = -0.4,
         y = cond_pow*100,
         li = cond_powCIl*100,
         ui = cond_powCIu*100,
         col=2, lwd=2, add=TRUE, err="y",
         barcol = "darkgreen")
  
  par(mar = c(0,4,1,0))
  
  x = seq(-0.4, 0, by = 0.0001)
  dens.thetao <- dnorm(x, mean = study$thetahato, sd = study$se_o)
  
  plot(x, dens.thetao,
       main = "", xlab= "", ylab = "", lwd =2,
       type = "l",
       axes = F,
       ylim = c(0, 100),
       xlim = c(-0.4, 0),
       col = "lightgrey")
  
  
  par(mar = c(4.1,0,0,0))
  
  plot(density(condPower_sample)$y/100, density(condPower_sample)$x,
       axes = F,
       main = "", xlab= "", ylab = "", lwd =2,
       type = "l",
       xlim = c(0, 0.5),
       ylim = c(0, 1.01),
       col = "lightgrey")
  
  
  
  
  ################################################################################
  ######################  Figure 4: Estimate - margin ############################
  ################################################################################
  
  
  ## Plotting of "shrinkage"
  par(mfrow = c(1, 1), 
      mar = c(5.1, 5, 4, 2), 
      las = 1, 
      pty = "s")
  library(scales)
  zalpha = p2z(0.025, "one.sided")
  
  plot(NA, 
       xlim = c(-1.55, 0.15), 
       ylim = c(-1.55, 0.15),
       xlab = expression(hat(theta)[RCT] - log(margin)), 
       ylab = expression(hat(theta)[RWE] - log(margin)))
  
  
  CIo_l = res$thetahato - log(res$margin) - zalpha*(res$se_o)
  CIo_u = res$thetahato - log(res$margin) + zalpha*(res$se_o)
  
  CIr_l <- res$thetahatr - log(res$margin) - zalpha*(res$se_r)
  CIr_u = res$thetahatr - log(res$margin) + zalpha*(res$se_r)
  
  segments(x0 = CIo_l, x1 = CIo_u, y0 = res$thetahatr - log(res$margin), 
           y1 = res$thetahatr - log(res$margin), 
           col = "lightgrey")
  
  segments(y0 = CIr_l, y1 = CIr_u, x0 = res$thetahato - log(res$margin), 
           x1 = res$thetahato - log(res$margin), 
           col = "lightgrey")
  
  
  abline(a = 0, b = 1, col = "red", lty = 2)
  points((res$thetahato - log(res$margin)),
         (res$thetahatr - log(res$margin)), 
         pch = 19, 
         col = alpha("black", 0.6))


  ################################################################################
  ######################  Figure A1: prop of success ############################
  ################################################################################
  
  
  po <- z2p(res$zo, alternative = "less")
  
  threshold <- seq(0.001, 0.1, by = 0.001)
  
  succ_rate <- av_cond_power <-  matrix(NA, nrow = length(threshold), 
                                        ncol = 3)
  
  
  for(i in 1:length(threshold)){
    succ_rate[i, 1] <- mean(res$TTR_p < threshold[i])
    succ_rate[i, 2] <- mean(res$golden_pS < threshold[i])
    succ_rate[i, 3] <- mean(res$controlled_pS < threshold[i])
    
    pow_2tr <- ifelse(po > threshold[i], 0, powerSignificance(zo = res$zo, c = res$se_o^2/res$se_r^2, level = threshold[i], 
                                                              alternative = "one.sided", designPrior = "conditional"))
    av_cond_power[i, 1] <- mean(pow_2tr)
    
    av_cond_power[i, 2] <- mean(powerReplicationSuccess(zo = res$zo, c = res$se_o^2/res$se_r^2, level = threshold[i], 
                                                        alternative = "one.sided", designPrior = "conditional", 
                                                        type = "golden"))
    
    av_cond_power[i, 3] <- mean(powerReplicationSuccess(zo = res$zo, c = res$se_o^2/res$se_r^2, level = threshold[i], 
                                                        alternative = "one.sided", designPrior = "conditional", 
                                                        type = "controlled"))
    
  }
  
  
  succrate25 <- succ_rate[threshold == 0.025]
  condpow25 <- av_cond_power[threshold == 0.025]
  
  par(las = 1, mfrow = c(1,1))
  matplot(threshold, succ_rate[, c(1,3)]*100, type = "l", lty = 1, 
          lwd = 2, 
          xlab = expression(alpha), 
          ylab = "%", 
          ylim = c(35, 100), col = c("black", "orange", 4)[c(1,3)])
  
  
  matlines(threshold, av_cond_power[, c(1, 3)]*100, lwd = 2, lty = 2, 
           col = c("black", "orange", 4)[c(1, 3)])
  
  
  segments(x0 = 0.025, y0 = succrate25[c(1,3)]*100, x1 = -0.5, 
           y1 = succrate25[c(1,3)]*100, lwd = 2, lty = 2, col = "lightgrey")
  
  segments(x0 = 0.025, y0 = 0, x1 = 0.025, 
           y1 = condpow25[c(1,3)]*100, lty = 2, lwd = 2, col = "lightgrey")
  
  segments(x0 = 0.025, y0 = condpow25[c(1,3)]*100, x1 = -0.5, 
           y1 = condpow25[c(1,3)]*100, lty = 2, lwd = 2, col = "lightgrey")
  
  
  axis(1, 0.025, 0.025, col = "red", col.axis = "red", cex.axis = 0.75)
  
  legend("bottomright", 
         c("Proportion of success", "Average conditional power"), 
         lty = c(1, 2), 
         lwd = 2, bty = "n", 
         cex = 0.9, 
         bg = "white")
  legend("bottomleft", 
         c("2TR", expression(paste("golden ", p[S])),
           expression(paste("controlled ", p[S])))[c(1,3)], 
         lty = 1, 
         col = c("black", "orange", 4)[c(1,3)], 
         lwd = 2,
         cex = 0.9, 
         bty = "n")
  
  
  ################################################################################
  ######################  Table 1 ############################
  ################################################################################
  
  
  df <- data.frame(Study = res$Study, round(res[,c(17, 18,  8, 7, 6)],5))
  df[,c(2, 3,  5, 6)] <- sapply(df[,c(2, 3, 5, 6)],biostatUZH::formatPval)
  colnames(df) <- c("Study", "$p_{RCT}$", "$p_{RWE}$","$c$", "$p_{TTR}$", "controlled $p_S$")

  df  
  
  
  
  ################################################################################
  ######################  Table 2 ############################
  ################################################################################
  
  # Calculation of the power
  
  power2TR_cond <- ifelse(z2p(res$zo, alternative = "less") < 0.025, 
                          powerSignificance(zo = res$zo, c = res$varianceRatio, 
                                            level = 0.025, designPrior = "conditional", 
                                            alternative = "one.sided"), 
                          0)
  
  power2TR_pred <- ifelse(z2p(res$zo, alternative = "less") < 0.025, 
                          powerSignificance(zo = res$zo, c = res$varianceRatio, 
                                            level = 0.025, designPrior = "predictive", 
                                            alternative = "one.sided"), 
                          0)
  
  powerRS_cond <- powerReplicationSuccess(zo = res$zo, c = res$varianceRatio, 
                                          level = 0.025, designPrior = "conditional", 
                                          alternative = "one.sided", 
                                          type = "controlled")
  
  powerRS_pred <- powerReplicationSuccess(zo = res$zo, c = res$varianceRatio, 
                                          level = 0.025, designPrior = "predictive", 
                                          alternative = "one.sided", 
                                          type = "controlled")

  
  # Make table with all types of power 
  
  power_mat <-  matrix(c(res$Study,
                         biostatUZH::formatPval(z2p(res$zo, alternative = "less")),
                         power2TR_cond*100, 
                         powerRS_cond*100, 
                         power2TR_pred*100,
                         powerRS_pred*100), ncol = 6)
  power_mat <- as.data.frame(power_mat)


  av_pow <- c("Average",  "",  mean(power2TR_cond*100), mean(powerRS_cond*100), 
              mean(power2TR_pred*100), mean(powerRS_pred*100))
  power_mat_ext <- rbind(power_mat, av_pow)
  
  # make numeric
  cols.num <- c("V3", "V4", "V5", "V6")
  power_mat_ext[cols.num] <- sapply(power_mat_ext[cols.num],as.numeric)
  
  
  colnames(power_mat_ext) <- c("Study", "$p_{RCT}$", 
                               "two-trials rule", "sceptical $p$", 
                               "two-trials rule", "sceptical $p$")
  
  
  power_mat_ext

  
