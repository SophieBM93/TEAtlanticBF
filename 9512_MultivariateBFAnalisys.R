##########################################################
# Downcore Benthic Foraminifera quantitative analyses    #
# for site GeoB9512-5                                    #
# Sofia Barragan-Montilla                                #
# Last modified: 11-1-2023                               # 
# Runned: 11-1-2023                                      # 
# For RStudio 2023.06.1 Build 524                        #
##########################################################

########################################################
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
setwd("C:/Users/Sofia Barragan Monti/Documents/PhDBFMultivariate/TEAtlanticBF")

# 3) Import Data for Diversity Calculations
ab1 <- read.csv("9512Sp.csv", 
                header = TRUE, row.names = 1, check.names = F)                     #Import data
ab1 <- replace(ab1,is.na(ab1),0)                                                   #Replacing NA values by 0
ab1<- ab1[rowSums(ab1[,2:148])>0,]                                                 #Delete 0 rows

######################
##Intervals (from Age Model)
YD <- tibble(ymin = 12.9, ymax =11.7, xmin = -Inf, xmax = Inf)
# HS1A <- tibble(ymin = 17.9, ymax =15.7, xmin = -Inf, xmax = Inf)
# HS1B <- tibble(ymin = 15.7, ymax =14.7, xmin = -Inf, xmax = Inf)
HS1 <- tibble(ymin = 18.36, ymax =15.47, xmin = -Inf, xmax = Inf)
HS2 <- tibble(ymin = 26.02, ymax =23.8, xmin = -Inf, xmax = Inf)
LGM <- tibble(ymin = 23, ymax =19, xmin = -Inf, xmax = Inf)
# CRobert<- tibble(ymin = 11.739, ymax =15.307, xmin = -Inf, xmax = Inf)

##################### 
### 4) Diversity Analyses Vegan Package

s <- specnumber(ab1[,])                  #Calculates simple Sp. Richness
H <- diversity(ab1[,], 
               index = "shannon", 
               MARGIN = 1, 
               base = exp(1))            #Calculates Shannon
A <- fisher.alpha(ab1[,])                #Calculates Fisher

Div <- data.frame(row.names(ab1),s,H,A)  #Creates df of calculated Div. values
colnames(Div)<- c("Age (ka)", "Richness", "Shannon Index", "Fisher Alpha Index") 
Div$`Age (ka)` <- as.numeric(Div$`Age (ka)`) #Make varibles numeric
Div$`Richness` <- as.numeric(Div$`Richness`) #Make varibles numeric
str(Div)                                     #Check types of variables
# write.csv(Div, file= "Diversity.csv")        #Saves calculations as csv

#Diveristy plots
dPlot <- ggplot(Div, aes(`Fisher Alpha Index`, `Age (ka)`))+
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
            data = LGM, alpha = 0.6, fill = "#c7c7c7", inherit.aes = FALSE)+
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
            data = HS1, alpha = 0.6, fill = "#c7c7c7",inherit.aes = FALSE)+
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
            data = HS2, alpha = 0.9, fill = "#c7c7c7",inherit.aes = FALSE)+
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
            data = YD, alpha = 0.6, fill = "#c7c7c7", inherit.aes = FALSE)+
  geom_lineh()+
  geom_lineh(aes(`Shannon Index`, `Age (ka)`), colour="#CC6666")+
  scale_y_reverse(guide="prism_minor",
                  breaks = seq(50,0,by=-5),
                  limits = c(27,0),
                  expand = c(0,0),
                  minor_breaks=seq(50, 0, -1))+
  scale_x_continuous(guide="prism_minor",
                  breaks = seq(0,20,by=5),
                  limits = c(0,20),
                  expand = c(0,0),
                  minor_breaks=seq(0,20, 1))+
  labs(title = "Benthic Foraminifera Diversity (Fisher Alpha")+
  ylab ("") +
  xlab ("") +
  scale_color_continuous(type = "viridis")+ 
  theme_minimal_vgrid()+
  theme(axis.ticks.y = element_blank(), 
        axis.line.y = element_blank(),
        axis.text.y = element_blank())
dPlot

######################
# 5) Calculates Relative Abundances and filter of most abundant species
#    later used for multivariate analyses

###Relative abundaces % 
ab1 <- ab1[rowSums(ab1[])>0,]                               #Deletes 0 rows
ab1P <- data.frame(ab1/rowSums(ab1)*100, check.names = F)   #Calculates %
rowSums(ab1P)                                               #Checks calculus went ok

### Filter species >2% in all samples
ab1P <- na.omit(ab1P)                                       #erases empty rows
crit <- c(count_col_if(gt(5), ab1P))                        #Criteria to filter (>2%)
crit2 <- c(ifelse(crit >= 1, 101,0))                        #Criteria to filter (in more than 1 sample)
ab1Pcr <- rbind(ab1P, crit2)                                #merge criteria with df
ab1PF <- select_if(ab1Pcr, (crit2==101))                    #Filter df
ab1PF <- ab1PF[-c(101), ]                                   #erase criteria
write.csv(ab1PF , file= "MatrizFiltrada%_5.csv")            #saves the resulting matrix
rowSums(ab1PF)                                              #Checks calculus went ok

###################################
# 6) Plot relative abundances of most important species 
r.df <- ab1PF
ydepth <- as.numeric(row.names(r.df))
Zones2 <- c(0, 11.7, 12.9, 14.6, 15.47, 18.36, 19, 23, 23.8, 26.02)

Stratiplot(r.df, ydepth,type = "poly", sort = "wa", rev=T, 
           varTypes="relative", zones = Zones2, 
           yticks = c(0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30),
           colour="black")
#Saves plot
pdf(file = "GeoB9512-5 BFAbund2.pdf", paper = "a4r", bg="transparent", height = 5, width = 30)
Stratiplot(r.df, ydepth,type = "poly", sort = "wa", rev=T, 
           varTypes="relative", zones = Zones2, 
           yticks = c(0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30),
           colour="black")
dev.off()

#####################
##  7) Multivariate analyses: Compute NMDS
ab1PF.Q <- (ab1PF*rowSums(ab1))/100                           #returns to BF counts
ab1PF.Q.Hel <-decostand(ab1PF.Q, method="hellinger")          #standarizes data

nmds <- metaMDS(ab1PF.Q.Hel, distance = "bray", k=7, 
                trymax=100000, autotransform = F)             #run NMDS
stressplot(nmds)                                              #check Stress Plot

# Import Environmental Data
env.9512 <- read.csv("9512EA.csv", check.names = F)
env.9512 <- replace(env.9512,is.na(env.9512),0)                                                  #Replacing NA values by 0
env.9512<- env.9512[rowSums(env.9512[,3:13])>0,]                                                 #Delete 0 rows

# Plot the NMDS 
fort <- fortify(nmds)
q1 <- ggplot()+
  geom_point(data=subset(fort, score=='sites'),
             mapping=aes(x=NMDS1, y=NMDS2, colour=factor(env.9512$Time)),
             alpha=0.5)+ 
  geom_segment(data=subset(fort, score=='species'), 
               mapping=aes(x=0, y=0, xend=NMDS1, yend=NMDS2),
               arrow = arrow(length =unit(0.015, "npc"), type="closed"), 
               colour="darkgray", size=0.5)+
  geom_text(data=subset(fort, score=='species'), 
            mapping = aes(label=label, x=NMDS1*1.1, y=NMDS2*1.1))+
  geom_abline(intercept = 0, slope=0, linetype="dashed", size=0.5, colour="black")+
  geom_vline(aes(xintercept=0), linetype="dashed", size=0.5, colour="black")+
  scale_x_continuous(guide="prism_minor",
                     breaks = seq(-20,20,by=1),
                     limits = c(-1.5,1.5),
                     expand = c(0,0),
                     minor_breaks=seq(-20,20, 0.5))+
  scale_y_continuous(guide="prism_minor",
                     breaks = seq(-20,20,by=0.5),
                     limits = c(-1,1),
                     expand = c(0,0),
                     minor_breaks=seq(-20,20, 0.25))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank(),
        axis.line = element_line(colour="black")) 
q1 + geom_text_repel() + geom_text() 

#### Make a 2 panel plot to reduce complexity
p1<- ggplot()+
  geom_point(data=subset(fort, score=='sites'),
             mapping=aes(x=NMDS1, y=NMDS2, colour=factor(env.9512$Time)),
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
                     breaks = seq(-20,20,by=1),
                     limits = c(-1.5,1.5),
                     expand = c(0,0),
                     minor_breaks=seq(-20,20, 0.5))+
  scale_y_continuous(guide="prism_minor",
                     breaks = seq(-20,20,by=0.5),
                     limits = c(-1,1),
                     expand = c(0,0),
                     minor_breaks=seq(-20,20, 0.25))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        legend.position = c(.25, 1),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(4,4,4,4),
        legend.title = element_blank(),
        legend.background = element_blank(),
        axis.line = element_line(colour="black")) +
  guides(color = guide_legend(ncol = 2))
p1

p2<- ggplot()+
  geom_point(data=subset(fort, score=='sites'),
             mapping=aes(x=NMDS1, y=NMDS2, colour=factor(env.9512$Time)),
             colour="black", alpha=0)+ 
  geom_segment(data=subset(fort, score=='species'), 
               mapping=aes(x=0, y=0, xend=NMDS1, yend=NMDS2),
               arrow = arrow(length =unit(0.015, "npc"), type="closed"), 
               colour="darkgray", size=0.5, alpha=0.5)+
  geom_text(data=subset(fort, score=='species'),
            mapping = aes(label=label, x=NMDS1*1.3, y=NMDS2*1.3),size=3, 
            family="serif", fontface="italic")+
  geom_abline(intercept = 0, slope=0, linetype="dashed", size=0.5, colour="black")+
  geom_vline(aes(xintercept=0), linetype="dashed", size=0.5, colour="black")+
  scale_x_continuous(guide="prism_minor",
                     breaks = seq(-20,20,by=1),
                     limits = c(-1.5,1.5),
                     expand = c(0,0),
                     minor_breaks=seq(-20,20, 0.5))+
  scale_y_continuous(guide="prism_minor",
                     breaks = seq(-20,20,by=0.5),
                     limits = c(-1,1),
                     expand = c(0,0),
                     minor_breaks=seq(-20,20, 0.25))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour="black"))
p2

#Visualize NMDS
ggarrange(p1, p2, ncol=1)                                          

# Saves NMDS
pdf(file = "GeoB9512-5 NMDS.pdf", paper = "a4r", bg="transparent", height = 5, width = 30)
ggarrange(p1, p2, ncol=1)
dev.off()

## Environmental fit variables in NMDS
env2 <- env.9512

#Extract the varibales of interest
envNMDS <- cbind(env2$`Age [ka]`, env2$`DO3[mL/L]`, env2$`Infaunals[%]`,
                      env2$`Epibenthic[%]`, env2$`Calcareous[%]`, env2$`Porcellaneous[%]`,
                      env2$`d13C (corrected) [o/oo]`, env2$`dw [permil]`)
envNMDS <- as.data.frame(envNMDS)
colnames(envNMDS)<- c("Age (ka)", "BWOx[mL/L]", "Infaunals [%]",
                      "Epifaunals [%]", "Calcareous [%]", 
                      "Porcellaneous [%]","d13C", "d18Osw") 
str(envNMDS)

#Plot variables in the NMDS
env <- envNMDS[2:8]
en = envfit(nmds, env, permutations = 999, na.rm = TRUE)                 #Do the environmnetal fit

data.scores = as.data.frame(scores(nmds, display="species"))                                 
en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)
en_coord_cat = as.data.frame(scores(en, "factors")) * ordiArrowMul(en)

#Plot the results
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
             mapping=aes(x=NMDS1, y=NMDS2, colour=factor(env.9512$Time)),
             alpha=0.5)+ 
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
                     limits = c(-2,2),
                     expand = c(0,0),
                     minor_breaks=seq(-20,20, 0.5))+
  scale_y_continuous(guide="prism_minor",
                     breaks = seq(-20,20,by=0.5),
                     limits = c(-1,1),
                     expand = c(0,0),
                     minor_breaks=seq(-20,20, 0.25))+
  xlab("")+ylab("")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        legend.position = c(1,1),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.title = element_blank(),
        legend.background = element_blank(),
        axis.line = element_line(colour = "black"))
p3

ggarrange(p3, p2, ncol = 1)

pdf(file = "NMDS9512_EnvFit.pdf", bg="transparent", height = 7, width = 7)
ggarrange(p3, p2, ncol=1)
dev.off()

###################################
## 8) Plot downcores Environmental variables
ea1 <- read.csv("9512EA.csv", check.names = F)
ea1 <- replace(ea1,is.na(ea1),0)                                                  #Replacing NA values by 0
ea1 <- ea1[rowSums(ea1[,3:5])>0,]                                                 #Delete 0 rows

#Stress Species %
SSPlot <- ggplot(ea1, aes(`Stress Species[%]`,`Age [ka]`))+  
  labs(title="Stress species") +
  scale_y_reverse(guide="prism_minor",
                  breaks = seq(50,0,by=-5),
                  limits = c(27,0),
                  expand = c(0,0),
                  minor_breaks=seq(50, 0, -1))+
  scale_x_continuous(guide="prism_minor",
                     breaks = seq(0,100,by=20),
                     limits = c(0,101),
                     expand = c(0,0),
                     minor_breaks=seq(0,100, 10))+
  labs(title="Stress Species (%)") +
  xlab("")+ 
  ylab("")+
  geom_areah()+
  theme_minimal_vgrid()+
  theme(panel.grid = element_blank(),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.y = element_blank(), 
        axis.line.y = element_blank(),
        axis.text.y = element_blank())
SSPlot

#TROX Model
Microhabitat <- data.frame(x2 = c(ea1$`Infaunals[%]`, ea1$`Epibenthic[%]`),  
                           y2 = ea1$`Age [ka]`,     
                           Morphogroups = c(rep("Infaunal", nrow(ea1)),
                                            rep("Epifaunal", nrow(ea1))))

MicrohabPlot <- ggplot(Microhabitat, aes(x2, y2)) +
  geom_areah(aes(fill = Morphogroups)) +
  scale_fill_grey(breaks = c('Infaunal', 'Epifaunal'))+
  scale_y_reverse(guide="prism_minor",
                  breaks = seq(50,0,by=-5),
                  limits = c(27,0),
                  expand = c(0,0),
                  minor_breaks=seq(50, 0, -1))+
  scale_x_continuous(guide="prism_minor",
                     breaks = seq(0,100,by=20),
                     limits = c(0,100),
                     expand = c(0,0),
                     minor_breaks=seq(0,100, 10))+
  labs(title="Morphogroups (%)") +
  xlab("")+ 
  ylab("")+ 
  geom_vline(xintercept = 33, linetype=2, size=0.5)+ 
  geom_vline(xintercept = 66, linetype=2, size=0.5) +
  theme_minimal_vgrid()+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.ticks.y = element_blank(), 
        axis.line.y = element_blank(),
        axis.text.y = element_blank())

MicrohabPlot                                    #grey are infunals

# Test type
Test <- data.frame(x2 = c(ea1$`Calcareous[%]`, ea1$`Porcellaneous[%]`, ea1$`Agglutinated[%]`),  
                   y2 = ea1$`Age [ka]`,     
                   Test_Composition = c(rep("Calcareous", nrow(ea1)),
                                        rep("Porcelaneous", nrow(ea1)), 
                                        rep("Agglutinated", nrow(ea1)))) # Reshape data frame

TestPlot <- ggplot(Test, aes(x2, y2)) +
  geom_areah(aes(fill = Test_Composition)) +
  scale_fill_grey(breaks = c('Calcareous', 'Porcelaneous', 'Agglutinated'))+
  fill_palette(palette = "grey")+
  scale_y_reverse(guide="prism_minor",
                  breaks = seq(50,0,by=-5),
                  limits = c(27,0),
                  expand = c(0,0),
                  minor_breaks=seq(50, 0, -1))+
  scale_x_continuous(guide="prism_minor",
                     breaks = seq(0,100,by=20),
                     limits = c(0,100),
                     expand = c(0,0),
                     minor_breaks=seq(0,100, 10))+
  labs(title="Test Composition (%)") +
  xlab("%")+ 
  ylab("")+ theme_bw()+ 
  theme_minimal_vgrid()+
  theme(legend.position = "none",
        panel.grid = element_blank())
TestPlot                                         #Black are porcelaneous, dark gray Calcareous

### Oxic preference %
Oxygen <- data.frame(x2 = c(ea1$`Oxic[%]`, ea1$`Suboxic[%]`, ea1$`Dysoxic[%]`),  
                     y2 = ea1$`Age [ka]`,     
                     OxygenGroup = c(rep("Oxic", nrow(ea1)),
                                     rep("Suboxic", nrow(ea1)), 
                                     rep("Dysoxic", nrow(ea1)))) # Reshape data frame

OxygenPlot <- ggplot(Oxygen, aes(x2, y2)) +
  geom_areah(aes(fill=OxygenGroup)) +
  fill_palette(palette = "grey")+
  scale_y_reverse(guide="prism_minor",
                  breaks = seq(50,0,by=-5),
                  limits = c(27,0),
                  expand = c(0,0),
                  minor_breaks=seq(50, 0, -1))+
  scale_x_continuous(guide="prism_minor",
                     breaks = seq(0,100,by=20),
                     limits = c(0,101),
                     expand = c(0,0),
                     minor_breaks=seq(0,100, 10))+
  labs(title="Oxygen Group(%)") +
  xlab("")+ 
  ylab("")+ theme_minimal_vgrid()+   
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.ticks.y = element_blank(), 
        axis.line.y = element_blank(),
        axis.text.y = element_blank())
OxygenPlot                                    #Black are suboxic, dark gray oxic, light geay dysoxic

## Dissolve Oxygen concentration from BFOI
DOPlot2 <- ggplot(ea1, aes(`DO3[mL/L]`, `Age [ka]`))+      #Plots the enhanced BFOI
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
            data = LGM, alpha = 0.6, fill = "#c7c7c7", inherit.aes = FALSE)+
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
            data = HS1, alpha = 0.6, fill = "#c7c7c7",inherit.aes = FALSE)+
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
            data = HS2, alpha = 0.9, fill = "#c7c7c7",inherit.aes = FALSE)+
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
            data = YD, alpha = 0.6, fill = "#c7c7c7", inherit.aes = FALSE)+
  geom_vline(xintercept = 0.1, linetype=2)+ 
  geom_vline(xintercept = 0.3, linetype=2) +
  geom_vline(xintercept = 1.5, linetype=2) +
  geom_vline(xintercept = 2, linetype=2)+
  scale_y_reverse(guide="prism_minor",
                  breaks = seq(50,0,by=-5),
                  limits = c(27,0),
                  expand = c(0,0),
                  minor_breaks=seq(50, 0, -1))+
  scale_x_continuous(guide="prism_minor",
                     breaks = seq(0,6,by=1),
                     limits = c(0,6),
                     expand = c(0,0),
                     minor_breaks=seq(0,6, 0.5))+
  geom_lineh()+
  labs(title="Dissolved Oxygen Concentration [mL/L]") +
  xlab("")+ 
  ylab("")+  
  theme_minimal_vgrid()+
  theme(panel.grid = element_blank(),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.y = element_blank(), 
        axis.line.y = element_blank(),
        axis.text.y = element_blank())
 
DOPlot2

#Plot all together
plot_grid(TestPlot, MicrohabPlot, SSPlot, 
          OxygenPlot, DOPlot2, dPlot, nrow = 1)

pdf(file = "GeoB9512-5 Oxyg.pdf", paper = "a4r", bg="transparent", height = 5, width = 30)
plot_grid(SSPlot, OxygenPlot, DOPlot2, nrow = 1)
dev.off()



