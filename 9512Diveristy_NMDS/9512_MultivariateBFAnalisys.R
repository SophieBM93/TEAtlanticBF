##########################################################
# Downcore Benthic Foraminifera quantitative analyses    #
# for site GeoB9512-5                                    #
# Sofia Barragan-Montilla                                #
# Last modified: 11-8-2023                               # 
# Runned: 11-8-2023                                      # 
# For RStudio 2023.06.1 Build 524                        #
##########################################################

########################################################
### Preapre Workspace

# 1) Clear Memory
rm(list=ls()) 

##################### 
##Packages

# install.packages('dplyr')       ##data manipulation 
library("dplyr")
# install.packages("tidyverse")   ##data manipulation etc.
library("tidyverse")
# install.packages("vegan")       ##For Diversity
library(vegan)
# install.packages("factoextra")  ##clustering algorithms & visualization
library(factoextra) 
# install.packages("FactoMineR")  ##data manipulation
library("FactoMineR") 
# install.packages("remotes")                        
# remotes::install_github("paleolimbot/tidypaleo")
# remotes::install_github("gavinsimpson/ggvegan")
library(tidypaleo)
library(ggvegan)
# install.packages("rioja")     ##clustering algorithms
library("rioja")  
library("expss")
library("cowplot")
library("ggpubr")
library("ggvegan")
# install.packages("analogue") 
library(analogue)
library("ggprism")
# install.packages("ggrepel")
library("ggrepel")
# install.packages("expss")
library("expss")

# 2) Set Directory
getwd()
setwd("C:/Users/Sofia Barragan Monti/Documents/PhDBFMultivariate/TEAtlanticBF/9512Diveristy_NMDS")

# 3) Import Data for Diversity Calculations
ab1 <- read.csv("9512Sp.csv", header = TRUE, check.names = F, row.names = 1)       #Import data
ab1 <- replace(ab1,is.na(ab1),0)                                                   #Replacing NA values by 0
ab1<- ab1[rowSums(ab1[,2:148])>0,]                                                 #Delete 0 rows
ab1 <- ab1[,2:148]


######################
## PART 1: Benthic Foraminifera Diversity and distribution

# 1.1) Defines Key Climatic Periods (As reported in the Manuscript)

YD <- tibble(ymin = 11.7, ymax =12.9, xmin = -Inf, xmax = Inf)
HS1 <- tibble(ymin = 18.5, ymax =14.79, xmin = -Inf, xmax = Inf)
HS2 <- tibble(ymin = 26.12, ymax =23.73, xmin = -Inf, xmax = Inf)
LGM <- tibble(ymin = 23, ymax =19, xmin = -Inf, xmax = Inf)


# 1.2) Calculates Diversity - Vegan Package

s <- specnumber(ab1[,])                       #Calculates simple Sp. Richness
H <- diversity(ab1[,], 
               index = "shannon", 
               MARGIN = 1, 
               base = exp(1))                 #Calculates Shannon
A <- fisher.alpha(ab1[,])                     #Calculates Fisher

Div <- data.frame(row.names(ab1),s,H,A)       #Creates df of calculated Div. values
colnames(Div)<- c("Age (ka)", "Richness", "Shannon Index", "Fisher Alpha Index") 
Div$`Age (ka)` <- as.numeric(Div$`Age (ka)`)  #Make varibles numeric
Div$`Richness` <- as.numeric(Div$`Richness`)  #Make varibles numeric
str(Div)                                      #Check types of variables
# write.csv(Div, file= "Diversity.csv")        #Saves calculations as csv

# 1.3) Plot Diversity
dPlot <- ggplot(Div, aes(`Fisher Alpha Index`, `Age (ka)`))+
  geom_hline(yintercept = 6.6, linetype = 2, linewidth = 0.5, alpha=0.5)+
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            data = LGM, alpha = 0.6, fill = "#c7c7c7", inherit.aes = FALSE)+
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            data = HS1, alpha = 0.4, fill = "#A6B578",inherit.aes = FALSE)+
  geom_hline(yintercept = 16.79, linetype=2, linewidth=0.5, alpha=0.5)+
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            data = HS2, alpha = 0.4, fill = "#A6B578",inherit.aes = FALSE)+
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            data = YD, alpha = 0.6, fill = "#c7c7c7", inherit.aes = FALSE)+
  geom_lineh()+
  scale_y_reverse(guide="prism_minor",
                  breaks = seq(50,0,by=-5),
                  limits = c(27,0),
                  expand = c(0,0),
                  minor_breaks=seq(50, 0, -1))+
  scale_x_continuous(guide="prism_minor",
                     breaks = seq(0,20,by=5),
                     limits = c(5,20),
                     expand = c(0,0),
                     minor_breaks=seq(0,20, 1))+
  labs(title = "Benthic Foraminifera Diversity (Fisher Alpha")+
  ylab ("") +
  xlab ("") +
  scale_color_continuous(type = "viridis")+ 
  theme_classic()+
  theme(axis.line.x = element_line(linetype = 1, size = 0.5),
        axis.ticks.x = element_line(linetype = 1, size = 0.5),
        plot.title = element_text(size = 11, family = "sans", hjust = 0.5, face="bold"),
        axis.text = element_text(size = 10, family = "sans"))
dPlot

# 1.4) Calculate Relative Abundances and filter of most abundant species
#    (later used for multivariate analyses - NMDS)

# Relative abundaces % 
ab1 <- read.csv("9512Gr2.CSV", check.names = F)
row.names(ab1) <-ab1$`Age (kyrs)`
ab1 <- ab1 [,2:82]
ab1 <- replace(ab1,is.na(ab1),0)
ab1 <- ab1[rowSums(ab1[,2:81])>0,]                          #Deletes 0 rows
ab1P <- data.frame(ab1/rowSums(ab1)*100, check.names = F)   #Calculates %
rowSums(ab1P)                                               #Checks calculus went ok

### Filter species >5% in all samples
ab1P <- na.omit(ab1P)                                       #erases empty rows
crit <- c(count_col_if(gt(7.5), ab1P))                        #Criteria to filter (>5%)
crit2 <- c(ifelse(crit >= 1, 101,0))                        #Criteria to filter (in more than 1 sample)
ab1Pcr <- rbind(ab1P, crit2)                                #merge criteria with df
ab1PF <- select_if(ab1Pcr, (crit2==101))                    #Filter df
ab1PF <- ab1PF[-c(101), ]                                   #erase criteria
write.csv(ab1PF , file= "MatrizFiltrada%_5.csv")            #saves the resulting matrix
rowSums(ab1PF)                                              #Checks calculus went ok

# 1.5) Plot relative abundances of most important species 

r.df <- ab1PF

ea1_b <- cbind(row.names(r.df), r.df[])

melted <- melt(ea1_b, id.vars = "row.names(r.df)")
str(melted)
names(melted)[1] <- "Age"
melted$Age <- as.numeric(melted$Age)

AbundSp <- ggplot(melted, aes(x=value,y= Age, fill=variable))+
  geom_hline(yintercept = 6.6, linetype = 2, linewidth = 0.5, alpha=0.7, colour="darkgreen")+
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            data = LGM, alpha = 0.6, fill = "#c7c7c7", inherit.aes = FALSE)+
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            data = HS1, alpha = 0.4, fill = "#A6B578",inherit.aes = FALSE)+
  geom_hline(yintercept = 16.79, linetype=2, linewidth=0.5, alpha=0)+
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            data = HS2, alpha = 0.4, fill = "#A6B578",inherit.aes = FALSE)+
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            data = YD, alpha = 0.6, fill = "#c7c7c7", inherit.aes = FALSE)+
  geom_hline(yintercept = 7.12, linetype = 2, linewidth = 0.5, alpha=0.5)+
  geom_hline(yintercept = 15.71, linetype = 2, linewidth = 0.5, alpha=0.5)+
  geom_hline(yintercept = 16.18, linetype = 2, linewidth = 0.5, alpha=0.5)+
  geom_hline(yintercept = 23.76, linetype = 1, linewidth = 0.5, alpha=0.5)+
  geom_hline(yintercept = 18.28, linetype = 1, linewidth = 0.5, alpha=0.5)+
  geom_hline(yintercept = 17.64, linetype = 1, linewidth = 0.5, alpha=0.5)+
  geom_hline(yintercept = 14.19, linetype = 1, linewidth = 0.5, alpha=0.5)+
  geom_hline(yintercept = 12.97, linetype = 1, linewidth = 0.5, alpha=0.5)+
  geom_hline(yintercept = 11.62, linetype = 1, linewidth = 0.5, alpha=0.5)+
  geom_hline(yintercept = 10.98, linetype = 1, linewidth = 0.5, alpha=0.5)+
  geom_areah()+
  scale_y_reverse(guide="prism_minor",
                  breaks = seq(50,0,by=-5),
                  limits = c(27,0),
                  expand = c(0,0),
                  minor_breaks=seq(50, 0, -1))+
  scale_x_continuous(guide="prism_minor",
                     breaks = seq(0,100,by=20),
                     limits = c(0,75),
                     expand = c(0,0),
                     minor_breaks=seq(0,100, 10))+
  scale_fill_manual(values = c("black","black","black","black", "black","black", 
                               "grey","grey", "black","grey","black", "#FFBABA",
                               "black","black","black","black","grey", "black",
                               "grey","#FFBABA", "#FFBABA", "#FFBABA", "#AFCAE6", "gray", 
                               "black","black","black"))+
  theme_classic()+
  theme(legend.position = "none",
        axis.line = element_line(linetype = 1, size = 0.5),
        axis.ticks = element_line(linetype = 1, size = 0.5),
        plot.title = element_text(size = 11, family = "sans", hjust = 0.5, face="bold"),
        axis.text = element_text(size = 10, family = "sans"))+
  facet_wrap(~ variable, nrow=1)
AbundSp 

#Saves plot
pdf(file = "GeoB9512-5_BFAbund.pdf", paper = "a4r", bg="transparent", height = 5, width = 30)
AbundSp 
dev.off()


######################
## PART 2: NMDS
##  2.1) Multivariate analyses: Compute NMDS

ab1PF2 <- ab1PF                                               #Data for NMDS (see 1.4)
ab1PF.Q <- (ab1PF2*rowSums(ab1))/100                          #returns to BF counts
 
ab1PF.Q.Hel <-decostand(ab1PF.Q, method="hellinger")          #standarizes data

nmds <- metaMDS(ab1PF.Q.Hel, distance = "bray", k=4, 
                trymax=1000000, autotransform = F)             #run NMDS
nmds
stressplot(nmds)                                              #check Stress Plot

nmds$points
scores(nmds)

df1 <- cbind(ab1PF.Q.Hel, envNMDS)
df1_b <- ab1PF.Q.Hel[,1:17]
ano <- anosim(ab1PF, df1$`BWOx[mL/L]`, distance = "bray", permutations = 9999)
ano


# 2.2 Import Environmental Data (Fot Plot)
env.9512 <- read.csv("9512EA.csv", check.names = F)
env.9512 <- replace(env.9512,is.na(env.9512),0)                         #Replacing NA values by 0
env.9512<- env.9512[rowSums(env.9512[,4:5])>0,]                         #Delete 0 rows

# Plot the NMDS 
fort <- fortify(nmds)
#### Make a 2 panel plot to reduce complexity
p1<- ggplot()+
  geom_point(data=subset(fort, score=='sites'),
             mapping=aes(x=NMDS1, y=NMDS2, shape = factor(env.9512$OM), 
                         colour=factor(env.9512$Time)),
             size=3, alpha=0.5)+ 
  geom_segment(data=subset(fort, score=='species'), 
               mapping=aes(x=0, y=0, xend=NMDS1, yend=NMDS2),
               arrow = arrow(length =unit(0.015, "npc"), type="closed"), 
               colour="darkgray", size=0, alpha=0)+
  geom_text(data=subset(fort, score=='species'),
            mapping = aes(label=label, x=NMDS1*1.3, y=NMDS2*1.3),size=3, 
            family="serif", fontface="italic", alpha=0)+
  geom_abline(intercept = 0, slope=0, linetype="dashed", size=0.5, colour="black")+
  geom_vline(aes(xintercept=0), linetype="dashed", size=0.5, colour="black")+
  scale_x_continuous(guide="prism_minor",
                     breaks = seq(-20,20,by=0.5),
                     limits = c(-1,1),
                     expand = c(0,0),
                     minor_breaks=seq(-20,20, 0.25))+
  scale_y_continuous(guide="prism_minor",
                     breaks = seq(-20,20,by=0.5),
                     limits = c(-1,1),
                     expand = c(0,0),
                     minor_breaks=seq(-20,20, 0.25))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        legend.position = c(0.25, 1),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(4,4,4,4),
        legend.title = element_blank(),
        legend.background = element_blank(),
        axis.line = element_line(linetype = 1, size = 0.5),
        axis.ticks = element_line(linetype = 1, size = 0.5),
        plot.title = element_text(size = 11, family = "sans", hjust = 0.5, face="bold"),
        axis.text = element_text(size = 10, family = "sans")) +
  guides(color = guide_legend(ncol = 2))
p1


# Saves NMDS
pdf(file = "GeoB9512_NMDSV3.pdf", paper = "a4r", bg="transparent", height = 5, width = 30)
ggarrange(p1)   
dev.off()

## Environmental fit variables in NMDS
## 2.3 Set Data
env2 <- env.9512                                           #Data for NMDS 

# 2.4 Extract the varibales of interest
envNMDS <- cbind(env2$`Age [ka]`, env2$`DO3[mL/L]`, env2$`Infaunals[%]`,
                 env2$`Epibenthic[%]`, env2$`Calcareous[%]`, env2$`Porcellaneous[%]`,
                 env2$FAI, env2$`(2) porcentage - C. robertsonianus`)
envNMDS <- as.data.frame(envNMDS)
colnames(envNMDS)<- c("Age (ka)", "BWOx[mL/L]", "Infaunals [%]",
                      "Epifaunals [%]", "Calcareous [%]", 
                      "Porcellaneous [%]", "FAI", "Elevated Epifauna") 
str(envNMDS)

# 2.5 Calculate environmental fit for variables
env <- envNMDS[2:8]
en = envfit(nmds, env, permutations = 999, na.rm = TRUE)                

data.scores = as.data.frame(scores(nmds, display="species"))                                 
en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)
en_coord_cat = as.data.frame(scores(en, "factors")) * ordiArrowMul(en)

# 2.6 Plot the results
p3<- ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2))+
  geom_segment(data = en_coord_cont, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
              size =0.5, alpha = 0.5, colour = "grey30") +
  # geom_point(data = en_coord_cat, aes(x = NMDS1, y = NMDS2), 
  #            shape = "diamond", size = 4, alpha = 0.6, colour = "navy") +
  # geom_text(data = en_coord_cat, aes(x = NMDS1, y = NMDS2+0.04), 
  #           label = row.names(en_coord_cat), colour = "navy", fontface = "bold", size=4) + 
  geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
             fontface = "bold", size=4, label = row.names(en_coord_cont)) + 
  geom_point(data=subset(fort, score=='sites'),
             mapping=aes(x=NMDS1, y=NMDS2, colour=factor(env.9512$Time), shape=factor(env.9512$OM)),
             alpha=0.5, size=4)+ 
  geom_segment(data=subset(fort, score=='species'), 
               mapping=aes(x=0, y=0, xend=NMDS1, yend=NMDS2),
               arrow = arrow(length =unit(0.015, "npc"), type="closed"), 
               colour="darkgray", size=0, alpha=0)+
  geom_text(data=subset(fort, score=='sites'),
            mapping = aes(label=label, x=NMDS1*1, y=NMDS2*1),size=3, 
            family="sans", fontface="italic", alpha=0)+
  geom_abline(intercept = 0, slope=0, linetype="dashed", size=0.5, colour="black")+
  geom_vline(aes(xintercept=0), linetype="dashed", size=0.5, colour="black")+
  scale_x_continuous(guide="prism_minor",
                     breaks = seq(-20,20,by=1),
                     limits = c(-1,1),
                     expand = c(0,0),
                     minor_breaks=seq(-20,20, 0.5))+
  scale_y_continuous(guide="prism_minor",
                     breaks = seq(-20,20,by=0.5),
                     limits = c(-1.5,1.5),
                     expand = c(0,0),
                     minor_breaks=seq(-20,20, 0.25))+
  xlab("")+ylab("")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        legend.position = c(0.1,1),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.title = element_blank(),
        legend.background = element_blank(),
        axis.line = element_line(colour = "black"))
p3







#####################
##  7) Multivariate analyses: Compute NMDS with taxonomic groups
# 1) Import data

ab2 <- read.csv("9512Gr2.csv", check.names = F, row.names = 1)
ab2 <- replace(ab2,is.na(ab2),0)                                                   #Replacing NA values by 0
ab2<- ab2[rowSums(ab2[,2:81])>0,]                                                 #Delete 0 rows
ab2 <- ab2[,2:81]

######################
# 2) Calculate Relative Abundances and filter of most abundant species

##Relative abundaces %
ab2 <- read.csv("9512EA_NMDS.CSV", check.names = F)
row.names(ab2) <-ab2$`Age [ka]`
ab2 <- ab2 [,2:13]
ab2 <- replace(ab2,is.na(ab2),0)
ab2 <- ab2[rowSums(ab2[,2:11])>0,]                          #Deletes 0 rows
ab2P <- data.frame(ab2/rowSums(ab2)*100, check.names = F)   #Calculates %
rowSums(ab2P)                                               #Checks calculus went ok

### Filter species >5% in all samples
ab2P <- na.omit(ab2P)                                       #erases empty rows
crit <- c(count_col_if(gt(5), ab2P))                        #Criteria to filter (>5%)
crit2 <- c(ifelse(crit >= 1, 101,0))                        #Criteria to filter (in more than 1 sample)
ab2Pcr <- rbind(ab2P, crit2)                                #merge criteria with df
ab2PF <- select_if(ab2Pcr, (crit2==101))                    #Filter df
ab2PF <- ab2PF[-c(101), ]                                   #erase criteria
write.csv(ab2PF , file= "MatrizFiltrada%_5(Groups).csv")    #saves the resulting matrix
rowSums(ab2PF)                                              #Checks calculus went ok

ab1PF.Q2 <- (ab2P*rowSums(ab1))/100                         #returns to BF counts

# 3) Standarize Data
ab1PF.Q.Hel2 <-decostand(ab1PF.Q2, method="hellinger")           #standarizes data

# 4) Run NMDS
nmds2 <- metaMDS(ab1PF.Q.Hel2, distance = "bray", k=8, 
                trymax=100000, autotransform = F)           #run NMDS
stressplot(nmds2)                                           #check Stress Plot

# 5) Plot Results
# Import Environmental Data For Plot
env.9512 <- read.csv("9512EA.csv", check.names = F)
env.9512 <- replace(env.9512,is.na(env.9512),0)             #Replacing NA values by 0
env.9512<- env.9512[rowSums(env.9512[,4:5])>0,]             #Delete 0 rows

# Plot the NMDS 

fort2 <- fortify(nmds2)                                     #Frotify NMDS for ggplot


p1_b<- ggplot()+
  geom_point(data=subset(fort2, score=='sites'),
             mapping=aes(x=NMDS1, y=NMDS2, shape = factor(env.9512$OM), 
                         colour=factor(env.9512$Time)),
             size=3, alpha=0.5)+ 
  geom_segment(data=subset(fort2, score=='species'), 
               mapping=aes(x=0, y=0, xend=NMDS1, yend=NMDS2),
               arrow = arrow(length =unit(0.015, "npc"), type="closed"), 
               colour="darkgray", size=0, alpha=0)+
  geom_text(data=subset(fort2, score=='species'),
            mapping = aes(label=label, x=NMDS1*1.3, y=NMDS2*1.3),size=3, 
            family="serif", fontface="italic", alpha=0)+
  geom_abline(intercept = 0, slope=0, linetype="dashed", size=0.5, colour="black")+
  geom_vline(aes(xintercept=0), linetype="dashed", size=0.5, colour="black")+
  scale_x_continuous(guide="prism_minor",
                     breaks = seq(-20,20,by=0.2),
                     limits = c(-1,1),
                     expand = c(0,0),
                     minor_breaks=seq(-20,20, 0.1))+
  scale_y_continuous(guide="prism_minor",
                     breaks = seq(-20,20,by=0.2),
                     limits = c(-1,1),
                     expand = c(0,0),
                     minor_breaks=seq(-20,20, 0.1))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        legend.position = c(0.3, 1),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(4,4,4,4),
        legend.title = element_blank(),
        legend.background = element_blank(),
        axis.line = element_line(linetype = 1, size = 0.5),
        axis.ticks = element_line(linetype = 1, size = 0.5),
        plot.title = element_text(size = 11, family = "sans", hjust = 0.5, face="bold"),
        axis.text = element_text(size = 10, family = "sans")) +
  guides(color = guide_legend(ncol = 2))
p1_b

# pdf(file = "GeoB9512_NMDS_GroupedData.pdf", paper = "a4r", bg="transparent", height = 5, width = 30)
# ggarrange(p1_b)   
# dev.off()

## Environmental fit variables in NMDS
## 2.3 Set Data
env2 <- env.9512                                           #Data for NMDS 

# 2.4 Extract the varibales of interest
envNMDS <- cbind(env2$`Age [ka]`, env2$`DO3[mL/L]`, env2$`Infaunals[%]`,
                 env2$`Epibenthic[%]`, env2$`Calcareous[%]`, env2$`Porcellaneous[%]`,
                 env2$FAI, env2$`(2) porcentage - C. robertsonianus`, env2$`Current Velocities (2)`)
envNMDS <- as.data.frame(envNMDS)
# envNMDS <- data.frame(lapply(envNMDS,as.numeric))
colnames(envNMDS)<- c("Age (ka)", "BWOx[mL/L]", "Infaunals [%]",
                      "Epifaunals [%]", "Calcareous [%]", 
                      "Porcellaneous [%]", "FAI", "Elevated Epifauna", "Current Velocities") 

# 2.5 Calculate environmental fit for variables
env <- envNMDS[2:9]
en2 = envfit(nmds2, env, permutations = 999, na.rm = TRUE)                

data.scores2 = as.data.frame(scores(nmds2, display="species"))                                 
en_coord_cont2 = as.data.frame(scores(en2, "vectors")) * ordiArrowMul(en2)
en_coord_cat2 = as.data.frame(scores(en2, "factors")) * ordiArrowMul(en2)

# 2.6 Plot the results
p3_b<- ggplot(data = data.scores2, aes(x = NMDS1, y = NMDS2))+
  geom_segment(data = en_coord_cont2, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               size =0.5, alpha = 0.5, colour = "grey30") +
  geom_text(data = en_coord_cont2, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
            fontface = "bold", size=4, label = row.names(en_coord_cont2)) + 
  geom_point(data=subset(fort2, score=='sites'),
             mapping=aes(x=NMDS1, y=NMDS2, colour=factor(env.9512$Time), shape=factor(env.9512$OM)),
             alpha=0.5, size=4)+ 
  geom_segment(data=subset(fort2, score=='species'), 
               mapping=aes(x=0, y=0, xend=NMDS1, yend=NMDS2),
               arrow = arrow(length =unit(0.015, "npc"), type="closed"), 
               colour="darkgray", size=0, alpha=0)+
  geom_text(data=subset(fort2, score=='sites'),
            mapping = aes(label=label, x=NMDS1*1, y=NMDS2*1),size=3, 
            family="sans", fontface="italic", alpha=0)+
  geom_abline(intercept = 0, slope=0, linetype="dashed", size=0.5, colour="black")+
  geom_vline(aes(xintercept=0), linetype="dashed", size=0.5, colour="black")+
  scale_x_continuous(guide="prism_minor",
                     breaks = seq(-20,20,by=1),
                     limits = c(-1,1),
                     expand = c(0,0),
                     minor_breaks=seq(-20,20, 0.5))+
  scale_y_continuous(guide="prism_minor",
                     breaks = seq(-20,20,by=0.5),
                     limits = c(-2,2),
                     expand = c(0,0),
                     minor_breaks=seq(-20,20, 0.25))+
  xlab("")+ylab("")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        legend.position = c(0.1,1),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.title = element_blank(),
        legend.background = element_blank(),
        axis.line = element_line(colour = "black"))
p3_b




# 1.3) Plot EBFOI
EBFOIPlot <- ggplot(env.9512, aes(EBFOI, `Age [ka]`))+
  geom_hline(yintercept = 6.6, linetype = 2, linewidth = 0.5, alpha=0.5)+
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            data = LGM, alpha = 0.6, fill = "#c7c7c7", inherit.aes = FALSE)+
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            data = HS1, alpha = 0.4, fill = "#A6B578",inherit.aes = FALSE)+
  geom_hline(yintercept = 16.79, linetype=2, linewidth=0.5, alpha=0.5)+
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            data = HS2, alpha = 0.4, fill = "#A6B578",inherit.aes = FALSE)+
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            data = YD, alpha = 0.6, fill = "#c7c7c7", inherit.aes = FALSE)+
  geom_lineh()+
  geom_vline(xintercept = 50, linetype=2, linewidth=0.5, alpha=0.5)+
  scale_y_reverse(guide="prism_minor",
                  breaks = seq(50,0,by=-5),
                  limits = c(27,0),
                  expand = c(0,0),
                  minor_breaks=seq(50, 0, -1))+
  scale_x_continuous(guide="prism_minor",
                     breaks = seq(0,100,by=20),
                     limits = c(0,80),
                     expand = c(0,0),
                     minor_breaks=seq(0,100, 10))+
  labs(title = "EBFOI")+
  ylab ("") +
  xlab ("") +
  scale_color_continuous(type = "viridis")+ 
  theme_classic()+
  theme(axis.line.x = element_line(linetype = 1, size = 0.5),
        axis.ticks.x = element_line(linetype = 1, size = 0.5),
        plot.title = element_text(size = 11, family = "sans", hjust = 0.5, face="bold"),
        axis.text = element_text(size = 10, family = "sans"))
EBFOIPlot

# 1.3) Plot EBFOI
BWO_Plot <- ggplot(env.9512, aes(`DO3[mL/L]`, `Age [ka]`))+
  geom_hline(yintercept = 6.6, linetype = 2, linewidth = 0.5, alpha=0.5)+
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            data = LGM, alpha = 0.6, fill = "#c7c7c7", inherit.aes = FALSE)+
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            data = HS1, alpha = 0.4, fill = "#A6B578",inherit.aes = FALSE)+
  geom_hline(yintercept = 16.79, linetype=2, linewidth=0.5, alpha=0.5)+
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            data = HS2, alpha = 0.4, fill = "#A6B578",inherit.aes = FALSE)+
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            data = YD, alpha = 0.6, fill = "#c7c7c7", inherit.aes = FALSE)+
  geom_lineh()+
  geom_point()+
  geom_vline(xintercept = 1.5, linetype=2, linewidth=0.5, alpha=0.5)+
  geom_vline(xintercept = 3, linetype=2, linewidth=0.5, alpha=0.5)+
  geom_vline(xintercept = 0.3, linetype=2, linewidth=0.5, alpha=0.5)+
  geom_vline(xintercept = 0.1, linetype=2, linewidth=0.5, alpha=0.5)+
  scale_y_reverse(guide="prism_minor",
                  breaks = seq(50,0,by=-5),
                  limits = c(27,0),
                  expand = c(0,0),
                  minor_breaks=seq(50, 0, -1))+
  scale_x_continuous(guide="prism_minor",
                     breaks = seq(0,100,by=1),
                     limits = c(0,6),
                     expand = c(0,0),
                     minor_breaks=seq(0,100, 0.5))+
  labs(title = "BWO")+
  ylab ("") +
  xlab ("") +
  scale_color_continuous(type = "viridis")+ 
  theme_classic()+
  theme(axis.line.x = element_line(linetype = 1, size = 0.5),
        axis.ticks.x = element_line(linetype = 1, size = 0.5),
        plot.title = element_text(size = 11, family = "sans", hjust = 0.5, face="bold"),
        axis.text = element_text(size = 10, family = "sans"))
BWO_Plot

ggarrange(EBFOIPlot, BWO_Plot, nrow=1)
