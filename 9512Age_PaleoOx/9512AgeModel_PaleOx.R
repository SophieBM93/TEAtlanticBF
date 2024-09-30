##########################################################
# Downcore Age model & Proxy uncertainty for GeoB9512-5  #
# Stefan Mulitza - Sofia Barragan-Montilla               #
# Last modified: 11-8-2023                               # 
# Runned: 11-8-2023                                      # 
# For R Version 1.4.1103                                 #
##########################################################

#####The age modelling code was adapted using the following references
####Blaauw, M., & J. A. Christen (2011). Flexible paleoclimate 
#age-depth models using an autoregressive gamma process. Bayesian Analysis, 6, 457-474.

###Heaton, T. J. et al. Marine20—The Marine Radiocarbon Age
#Calibration Curve (0–55,000 cal BP). Radiocarbon 62, 779–820 (2020).

########################################################
# Part 1: Age Modelling

# 1) Clear Memory
rm(list=ls()) 

# 2) Set working Directory (costumize for your own path)
setwd("C:/Users/Sofia Barragan Monti/Documents/PhDBFMultivariate/TEAtlanticBF/9512Age_PaleoOx") 

# 3) Load Packages
library(rbacon) #Load BACON & IntCal libraries
library(IntCal)
library(ggplot2)
library(tidyverse)
# install.packages("remotes")                        
# remotes::install_github("paleolimbot/tidypaleo")
library(tidypaleo)
library(cowplot)
library(datasets)

########################################################
# 4) Set core limits
d.min=2.5      #Minimum depth
d.max=547.5    #Maximum depth
nmc=10000      #Number of iterations (age models and time series to generate)

########################################################
# 5) Load Proxy data
#14C Data (Age data is saved in a separate folder called Bacon_runs)
#Radiocarbon dates from https://doi.org/10.1594/PANGAEA.962899 
Age<-read.csv(file="C:/Users/Sofia Barragan Monti/Documents/PhDBFMultivariate/TEAtlanticBF/9512Age_Proxies/Bacon_runs/GeoB9512-5/GeoB9512-5.csv", 
                  header=TRUE, stringsAsFactors=FALSE,sep=",",na.strings="-", check.names = F)

#Paleo-Oxyegenation Record GeoB9512-5 (saved in this repository) 
BWOx<-read.csv(file="9512_EA.csv", 
              header=TRUE, stringsAsFactors=FALSE,sep=",",na.strings="-", check.names = F)
BWOx<-data.frame(BWOx$`Depth [cm]`, BWOx$`Inferred Dissolve Oxygen Concentration [umol/kg]`)
colnames(BWOx)<-c("Depth [cm]", "BWOx")
BWOx <- replace(BWOx,is.na(BWOx),0)                                                         #Replacing NA values by 0
BWOx<- BWOx[(BWOx[,2])>0,]                                                                  #Delete 0 rows

########################################################
# 6) Run age model & save age model data (Figure 4)

#Reload old run
Bacon("GeoB9512-5",ask=FALSE,ssize=1.5*nmc,d.min=d.min,d.max=d.max,acc.mean=50, 
      mem.mean=0.7,mem.strength=5,acc.shape=1.5, run=FALSE)

agedepth()                                                                                  #loads previous model
# Run Bacon (Run this lines (65-68) to re-calculate the Age model)
# Bacon(core="GeoB9512-5", ask=FALSE,ssize=1.5*nmc,d.min=d.min,
#       d.max=d.max,acc.mean=50,d.by = 1, mem.mean=0.7,mem.strength=5,
#       acc.shape=1.5,t.a=20,t.b=21, res=20)

# Extract ages cm resolution
depthinCm=matrix(seq(d.min,d.max, 1))
AllAge<-sapply(depthinCm, Bacon.Age.d) # Extract ages with Bacon
NoCol<-dim(AllAge)[1] # Find number of age models

# Bacon produces more age models than requested, we only keep the defined number (nmc)
AllChron<-t(AllAge)[,(NoCol-nmc+1):NoCol] # Restrict age models to last nmc 
AllChron.median<-matrix(apply(AllChron,1,median)) # Calculate median age
AllChron.mean<-matrix(apply(AllChron,1,mean)) # Calculate median age

# Calculate 95% and 68% age envelope & save age model
AllChron.Qt<-matrix(0,nrow=length(depthinCm[,1]),ncol=4,byrow=FALSE)
probs<-c(0.025,0.16,0.84,0.975)
AllChron.Qt<-t(apply(AllChron,1,quantile,probs=probs,na.rm=TRUE))
dfAllChron<-data.frame(cbind(depthinCm,AllChron.median,AllChron.Qt)) # Merge fields
colnames(dfAllChron)<-c("depth_cm","Median Age","Q0.025","Q0.16", "Q0.84", "Q0.975")
write.csv(dfAllChron,file=sprintf("OPGeoB9512-5_AgeModel2.csv"), row.names=FALSE)


########################################################
# 7) Write out all simulated ages for the proxy depths
#BWOx
BWOxMCAG<-sapply(BWOx[,1], Bacon.Age.d) # Extract ages with Bacon for Isotopes
BWOxneff<-dim(BWOxMCAG)[1] # Determine number of columns
BWOxchron.dat<-t(BWOxMCAG)[,(BWOxneff-nmc+1):BWOxneff] # Tailor to desired number of iterations (nmc)
BWOxchron.median<-apply(BWOxchron.dat,1,median) # Calculate median age

########################################################
# 8) Add uncertainty to proxy
shapiro.test(BWOx$BWOx)                                                    # Normality test


BWOxdat.unc<-matrix(0,nrow=length(BWOx[,2]),ncol=nmc,byrow=FALSE) # Prepare matrix for d18O + noise
BWOxdat.prx<-matrix(BWOx[,2],nrow=length(BWOx[,2]),ncol=nmc,byrow=FALSE) # Prepare matrix with data
BWOxsig.prx<-matrix(0,nrow=length(BWOx[,2]),ncol=nmc,byrow=FALSE) # Prepare matrix for noise
for(j in 1:length(BWOx[,2])){BWOxsig.prx[j,]<-rnorm(n=nmc,sd=7.07)} # Put noise into matrix

BWOxdat.unc<-BWOxdat.prx+BWOxsig.prx # Add noise to data and store in matrix

BWOxdat.prx<-matrix(BWOx[,2],nrow=length(BWOx[,2]),ncol=nmc,byrow=FALSE)   # Prepare matrix with data
BWOxdat.unc<-BWOxdat.prx                                                   # Add noise to data and store in matrix

########################################################
# 9) Interpolate BWOx ensembles to the median age

BWOxMedAge<-matrix(NA, nrow=length(BWOxchron.median), ncol=nmc)            # Matrix for interpolated values
for (i in 1:nmc)
{
BWOxMedAge[,i]<-approx(BWOxchron.dat[,i],BWOxdat.unc[,i],BWOxchron.median)$y 
}
BWOxMedAge.mean<-matrix(apply(BWOxMedAge,1,mean,na.rm=TRUE))               # Calculate mean BWOx

########################################################
# 10) Interpolate BWOx Ensemble to the median age

BWOxMedAgeHigh<-matrix(NA, nrow=length(BWOxchron.median), ncol=nmc) # Matrix for interpolated values
for (i in 1:nmc)
{
  BWOxMedAgeHigh[,i]<-approx(BWOxchron.dat[,i],BWOxdat.unc[,i],BWOxchron.median)$y 
}

BWOxMedAgeHigh.mean<-matrix(apply(BWOxMedAgeHigh,1,mean,na.rm=TRUE)) # Calculate mean d18O

########################################################
# 12) Extract quantiles from BWOx time series
probs <- c(0.025,0.16,0.84,0.975)
QtsBWOx <- t(apply(BWOxMedAgeHigh,1,quantile,probs=probs,na.rm=TRUE))

########################################################
# 13) Save data -Optional
dfBWOx<-data.frame(cbind(BWOx[,1],BWOxchron.median/1000,BWOxMedAgeHigh.mean,QtsBWOx)) # Merge fields
colnames(dfBWOx)<-c("Depth [cm]","Median Age [kyrs]","Mean BWOx [umol/kg]","Q0.025","Q0.16", "Q0.84", "Q0.975")
write.csv(dfBWOx,file=sprintf("OPGeoB9512-5_BWOxQts.csv"), row.names=FALSE)

################################
# 14) XRF to CORROBORATE Heinrich stadials and Younger Dryas
# data from https://doi.org/10.1029/2018PA003359
XRF<-read.csv(file="GeoB9512-5_XRF.csv",
              header=TRUE, stringsAsFactors=FALSE,sep=",",na.strings="-", check.names = F)

# Write out all simulated ages for depths
XRFMCAG<-sapply(XRF[,1], Bacon.Age.d) # Extract ages with Bacon
XRFneff<-dim(XRFMCAG)[1] # Determine number of columns
XRFchron.dat<-t(XRFMCAG)[,(XRFneff-nmc+1):XRFneff] # Tailor to desired number of iterations (nmc)
XRFchron.median<-apply(XRFchron.dat,1,median) # Calculate median age

XRF_Age<-data.frame(cbind((XRFchron.median/1000),XRF)) # Merge fields
colnames(XRF_Age)<-c("Median Age [kyrs]","Depth (cm)", "Ti/Ca",
                     "Fe/Ca", "Fe/K", "Mn", "Total Counts")
write.csv(XRF_Age,file=sprintf("OPGeoB9512-5_XRF_Age.csv"), row.names=FALSE)

########################################################
# Part 2: Results Plotting

########################################################
# 1) Load Packages
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
# install.packages("ggprism")
library("ggprism")
## Adding Graphic elements (Time Intervals)  **Callib with Fe/Ca

########################################################
# *If needed: load data

dfBWOx<-read.csv("OPGeoB9512-5_BWOxQts.csv", check.names = F)
XRF<-read.csv("OPGeoB9512-5_XRF_Age.csv", check.names=F)

ea1 <- read.csv("9512_EA.csv", check.names=FALSE)   #GeoB9512-5 Paleoenvironmental conditions from Benthic Foraminifera
ea1 <- replace(ea1,is.na(ea1),0)                    #Replacing NA values by 0
ea1 <- ea1[rowSums(ea1[,26:27])>0,]                 #Delete 0 rows

########################################################
# 1) Define Key Climatic Periods (As reported in the Manuscript)

YD <- tibble(ymin = 11.7, ymax =12.9, xmin = -Inf, xmax = Inf)
HS1 <- tibble(ymin = 18.5, ymax =14.79, xmin = -Inf, xmax = Inf)
HS2 <- tibble(ymin = 26.12, ymax =23.73, xmin = -Inf, xmax = Inf)
LGM <- tibble(ymin = 23, ymax =19, xmin = -Inf, xmax = Inf)

########################################################
# 2) Plot Elemental ratios (uncalibrated) with the new age model to confirm HS2 and HS1 
XRFplot <- ggplot(XRF, aes(`Fe/Ca`, `Median Age [kyrs]`))+
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            data = LGM, alpha = 0.6, fill = "#c7c7c7", inherit.aes = FALSE)+
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            data = HS1, alpha = 0.4, fill = "#A6B578",inherit.aes = FALSE)+
  geom_hline(yintercept = 16.79, linetype=2, linewidth=0.5, alpha=0.5)+
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            data = HS2, alpha = 0.4, fill = "#A6B578",inherit.aes = FALSE)+
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            data = YD, alpha = 0.6, fill = "#c7c7c7", inherit.aes = FALSE)+
  geom_lineh(XRF, mapping=aes(`Fe/Ca`, `Median Age [kyrs]`),
             colour="black", linetype = 1, linewidth = 0.5)+
  scale_y_reverse(guide="prism_minor",
                  breaks = seq(50,0,by=-5),
                  limits = c(27,0),
                  expand = c(0,0),
                  minor_breaks=seq(50, 0, -1))+
  scale_x_continuous(guide="prism_minor",
                     breaks = seq(-2,20,by=2),
                     limits = c(0,10),
                     expand = c(0,0),
                     minor_breaks=seq(-2,20, 1))+
  labs(title="Fe/Ca")+ xlab("")+ ylab("Age (Kyrs)")+
  theme_classic()+
  theme(axis.line.x = element_line(linetype = 1, linewidth = 0.5),
        axis.ticks.x = element_line(linetype = 1, linewidth = 0.5),
        plot.title = element_text(size = 11, family = "sans", hjust = 0.5, face="bold"),
        axis.text.x = element_text(size = 10, family = "sans"))
XRFplot

############################
# 3) Plot Figure 2 

##Ph/Th from the Atlantic Ocean https://www.nature.com/articles/nature02494
PhTh_A <- read.csv(file="PhTh_McManus.csv", check.names = F, header=TRUE, 
                   stringsAsFactors=FALSE,sep=",",na.strings="-")
colnames(PhTh_A)<-c("Depth [cm]", "Age [kyrs]", "Pa/Th238", "Error238", "Pa/Th232","Error232")

#https://www.nature.com/articles/nature14059
PhTh_C <- read.csv(file="Bohm2015.csv", check.names = F, header=TRUE, 
                   stringsAsFactors=FALSE,sep=",",na.strings="-")
colnames(PhTh_C)<-c("Event", "Depth [cm]", "Age [kyrs]", "231Pa/230Th", "2sd %", "2SD")

#Plot Atlantic Ph/Th
PhThPlot <- ggplot()+
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
  geom_lineh(PhTh_C, mapping=aes(`231Pa/230Th`, `Age [kyrs]`),
             colour="#CC6666", linetype = 1, size = 0.25)+
  geom_errorbar(PhTh_C,mapping=aes(y= `Age [kyrs]`,xmin = `231Pa/230Th` - `2SD`,
                                   xmax = `231Pa/230Th` + `2SD`),
                colour="#CC6666", linetype = 1, width = 0.1, size=0.1)+
  geom_errorbar(PhTh_A,mapping=aes(y= `Age [kyrs]`,xmin = `Pa/Th232` - Error232,
                                   xmax = `Pa/Th232` + Error232),
                colour= "#66CC99", linetype = 1, width = 0.5)+
  geom_lineh(PhTh_A, mapping=aes(`Pa/Th232`, `Age [kyrs]`),
             colour="forestgreen", linetype = 1, size = 0.25)+
  geom_errorbar(PhTh_A,mapping=aes(y= `Age [kyrs]`,xmin = `Pa/Th238` - Error238,
                                   xmax = `Pa/Th238` + Error238),
                colour= "darkgray", linetype = 1, width = 0.5)+
  geom_lineh(PhTh_A, mapping=aes(`Pa/Th238`, `Age [kyrs]`),
             colour="black", linetype = 1, size = 0.25)+
  scale_y_reverse(guide="prism_minor",
                  breaks = seq(180,0,by=-5),
                  limits = c(27,0),
                  expand = c(0,0),
                  minor_breaks=seq(180, 0, -1))+
  scale_x_reverse(guide="prism_minor",
                  breaks = seq(0.14,0,by=-0.02),
                  limits = c(0.11,0.05),
                  expand = c(0,0),
                  minor_breaks=seq(0.14,0,-0.01))+
  labs(title="231Pa/230Th")+
  xlab("")+ ylab("")+ theme_classic()+
  theme_classic()+
  theme(axis.line.x = element_line(linetype = 1, size = 0.5),
        axis.ticks.x = element_line(linetype = 1, size = 0.5),
        plot.title = element_text(size = 11, family = "sans", hjust = 0.5, face="bold"),
        axis.text.x = element_text(size = 10, family = "sans"))
PhThPlot

#SST Gradients 
##Load Data first
SST_Grad <- read.csv("SST_Gradients.csv", check.names = F)

#Plot SSTs Gradients
SSTPlot <- ggplot(SST_Grad, aes(Diff, `Age [ka BP] SU81`))+      #Plots the enhanced BFOI
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            data = LGM, alpha = 0.6, fill = "#c7c7c7", inherit.aes = FALSE)+
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            data = HS1, alpha = 0.4, fill = "#A6B578",inherit.aes = FALSE)+
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            data = HS2, alpha = 0.4, fill = "#A6B578",inherit.aes = FALSE)+
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            data = YD, alpha = 0.6, fill = "#c7c7c7", inherit.aes = FALSE)+
  scale_y_reverse(guide="prism_minor",
                  breaks = seq(50,0,by=-5),
                  limits = c(27,0),
                  expand = c(0,0),
                  minor_breaks=seq(50, 0, -1))+
  scale_x_continuous(guide="prism_minor",
                     breaks = seq(0,20,by=5),
                     limits = c(8,19),
                     expand = c(0,0),
                     minor_breaks=seq(0,20, 1))+
  geom_lineh()+
  geom_point(color="black", size=1, alpha=0.5)+
  labs(title="North Atlantic SST Gradient (°C)") +
  xlab("")+ 
  ylab("")+  
  theme_classic()+
  theme(axis.ticks.y = element_blank(), 
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.x = element_line(linetype = 1, size = 0.5),
        axis.ticks.x = element_line(linetype = 1, size = 0.5),
        plot.title = element_text(size = 11, family = "sans", hjust = 0.5, face="bold"),
        axis.text.x = element_text(size = 10, family = "sans"))
SSTPlot

#Plot Benthic Foraminifera Oxygen Pref 
##Reshape dataframe to plot
OxicPref <- data.frame(x2 = c(ea1$`Oxic[%]`, ea1$`Suboxic[%]`, ea1$`Dysoxic[%]`),  
                       y2 = ea1$`Age [ka]`,     
                       Test_Composition = c(rep("Oxic", nrow(ea1)),
                                            rep("Suboxic", nrow(ea1)), 
                                            rep("Dysoxic", nrow(ea1)))) # Reshape data frame
##Plot
OxicPrefPlot <- ggplot(OxicPref, aes(x2, y2)) +
  geom_areah(aes(fill = factor(Test_Composition, levels=c('Dysoxic', 'Suboxic', 'Oxic')))) +
  # scale_fill_grey(breaks = c('Oxic', 'Suboxic', 'Dysoxic'))+
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
  labs(title="Oxic preferences (%)") +
  xlab("")+ 
  ylab("")+ theme_bw()+ 
  theme_classic()+
  theme(legend.position = "none", 
        axis.ticks.y = element_blank(), 
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.x = element_line(linetype = 1, size = 0.5),
        axis.ticks.x = element_line(linetype = 1, size = 0.5),
        plot.title = element_text(size = 11, family = "sans", hjust = 0.5, face="bold"),
        axis.text = element_text(size = 10, family = "sans"))
OxicPrefPlot  

## Dissolve Oxygen concentration from EBFOI
BWOxPlot <- ggplot(dfBWOx, aes(`Mean BWOx [umol/kg]`, `Median Age [kyrs]`))+      #Plots the enhanced BFOI
  geom_hline(yintercept = 3.4, linetype = 2, linewidth = 0.5, alpha=0.5)+
  geom_hline(yintercept = 5.1, linetype = 2, linewidth = 0.5, alpha=0.5)+
  geom_hline(yintercept = 10.85, linetype = 2, linewidth = 0.5, alpha=0.5)+
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            data = LGM, alpha = 0.6, fill = "#c7c7c7", inherit.aes = FALSE)+
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            data = HS1, alpha = 0.4, fill = "#A6B578",inherit.aes = FALSE)+
  geom_hline(yintercept = 15.75, linetype=2, linewidth=0.5, alpha=0.5)+
  geom_hline(yintercept = 16.4, linetype=2, linewidth=0.5, alpha=0.5)+
  geom_hline(yintercept = 17.55, linetype=2, linewidth=0.5, alpha=0.5)+
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            data = HS2, alpha = 0.4, fill = "#A6B578",inherit.aes = FALSE)+
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            data = YD, alpha = 0.6, fill = "#c7c7c7", inherit.aes = FALSE)+
  scale_y_reverse(guide="prism_minor",
                  breaks = seq(50,0,by=-5),
                  limits = c(27,0),
                  expand = c(0,0),
                  minor_breaks=seq(50, 0, -1))+
  scale_x_continuous(guide="prism_minor",
                     breaks = seq(0,300,by=50),
                     limits = c(80,260),
                     expand = c(0,0),
                     minor_breaks=seq(0,300, 10))+
  geom_lineh(data=ea1, mapping=aes(`Inferred Dissolve Oxygen Concentration [umol/kg]`, `Age [ka]`),
             color="red", size=0.5, alpha=0.5)+
  geom_point(data=ea1, mapping=aes(`Inferred Dissolve Oxygen Concentration [umol/kg]`, `Age [ka]`), 
             color="red", size=1, alpha=0.5)+
  # geom_ribbon(dfBWOx,mapping=aes(xmin = Q0.16, xmax = Q0.84), fill = "black", alpha=0.2)+
  geom_ribbon(dfBWOx,mapping=aes(xmin = Q0.025, xmax = Q0.975), fill = "black", alpha=0.2)+
  geom_lineh(dfBWOx, mapping=aes(`Mean BWOx [umol/kg]`, `Median Age [kyrs]`), 
             colour="black", linetype = 1, size = 0.6)+
  labs(title="BWOx [mL/L]") +
  xlab("")+ 
  ylab("")+  
  theme_classic()+
  theme(axis.ticks.y = element_blank(), 
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.x = element_line(linetype = 1, size = 0.5),
        axis.ticks.x = element_line(linetype = 1, size = 0.5),
        plot.title = element_text(size = 11, family = "sans", hjust = 0.5, face="bold"),
        axis.text.x = element_text(size = 10, family = "sans"))
BWOxPlot

## Figure 2
ggarrange(PhThPlot, SSTPlot, OxicPrefPlot, BWOxPlot, nrow = 1)


############################
##
BWOxPlot_Fig3 <- ggplot(dfBWOx, aes(`Mean BWOx [umol/kg]`, `Median Age [kyrs]`))+      #Plots the enhanced BFOI
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            data = LGM, alpha = 0.6, fill = "#c7c7c7", inherit.aes = FALSE)+
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            data = HS1, alpha = 0.4, fill = "#A6B578",inherit.aes = FALSE)+
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            data = HS2, alpha = 0.4, fill = "#A6B578",inherit.aes = FALSE)+
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            data = YD, alpha = 0.6, fill = "#c7c7c7", inherit.aes = FALSE)+
  geom_hline(yintercept = 3.4, linetype = 2, linewidth = 0.5, alpha=0.5)+
  geom_hline(yintercept = 5.1, linetype = 2, linewidth = 0.5, alpha=0.5)+
  geom_hline(yintercept = 10.85, linetype = 2, linewidth = 0.5, alpha=0.5)+
  geom_hline(yintercept = 15.75, linetype=2, linewidth=0.5, alpha=0.5)+
  geom_hline(yintercept = 16.4, linetype=2, linewidth=0.5, alpha=0.5)+
  geom_hline(yintercept = 17.55, linetype=2, linewidth=0.5, alpha=0.5)+
  scale_y_reverse(guide="prism_minor",
                  breaks = seq(50,0,by=-5),
                  limits = c(20,10),
                  expand = c(0,0),
                  minor_breaks=seq(50, 0, -1))+
  scale_x_continuous(guide="prism_minor",
                     breaks = seq(0,300,by=50),
                     limits = c(80,260),
                     expand = c(0,0),
                     minor_breaks=seq(0,300, 10))+
  geom_lineh(data=ea1, mapping=aes(`Inferred Dissolve Oxygen Concentration [umol/kg]`, `Age [ka]`),
             color="red", size=0.5, alpha=0.5)+
  geom_point(data=ea1, mapping=aes(`Inferred Dissolve Oxygen Concentration [umol/kg]`, `Age [ka]`), 
             color="red", size=1, alpha=0.5)+
  geom_ribbon(dfBWOx,mapping=aes(xmin = Q0.025, xmax = Q0.975), fill = "black", alpha=0.2)+
  geom_lineh(dfBWOx, mapping=aes(`Mean BWOx [umol/kg]`, `Median Age [kyrs]`), 
             colour="black", linetype = 1, size = 0.6)+
  labs(title="BWOx [mL/L]") +
  xlab("")+ 
  ylab("")+  
  theme_classic()+
  theme(axis.ticks.y = element_blank(), 
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.x = element_line(linetype = 1, size = 0.5),
        axis.ticks.x = element_line(linetype = 1, size = 0.5),
        plot.title = element_text(size = 11, family = "sans", hjust = 0.5, face="bold"),
        axis.text.x = element_text(size = 10, family = "sans"))
BWOxPlot_Fig3

## 4) Figure 5

#Stress Species
SSPlot <- ggplot(ea1, aes(`Stress Species[%]`,`Age [ka]`))+
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
  labs(title="Stress species") +
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
  labs(title="Stress Species (%)") +
  xlab("")+ 
  ylab("")+
  geom_areah()+
  theme_classic()+
  theme(axis.ticks.y = element_blank(), 
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.x = element_line(linetype = 1, size = 0.5),
        axis.ticks.x = element_line(linetype = 1, size = 0.5),
        plot.title = element_text(size = 11, family = "sans", hjust = 0.5, face="bold"),
        axis.text = element_text(size = 10, family = "sans"))
SSPlot

Inf_Fig3 <- ggplot(ea1, aes(`Infaunals[%]`, `Age [ka]`))+      #Plots the enhanced BFOI
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            data = LGM, alpha = 0.6, fill = "#c7c7c7", inherit.aes = FALSE)+
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            data = HS1, alpha = 0.4, fill = "#A6B578",inherit.aes = FALSE)+
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            data = HS2, alpha = 0.4, fill = "#A6B578",inherit.aes = FALSE)+
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            data = YD, alpha = 0.6, fill = "#c7c7c7", inherit.aes = FALSE)+
  scale_y_reverse(guide="prism_minor",
                  breaks = seq(50,0,by=-5),
                  limits = c(27,0),
                  expand = c(0,0),
                  minor_breaks=seq(50, 0, -1))+
  scale_x_continuous(guide="prism_minor",
                     breaks = seq(0,100,by=20),
                     limits = c(30,100),
                     expand = c(0,0),
                     minor_breaks=seq(0,100, 10))+
  geom_lineh(data=ea1, mapping=aes(`Infaunals[%]`, `Age [ka]`),
             color="black", size=0.5, alpha=0.5)+
  geom_point(data=ea1, mapping=aes(`Infaunals[%]`, `Age [ka]`), 
             color="black", size=1, alpha=0.5)+
  labs(title="Infaunal Foraminifera [%]") +
  xlab("")+ 
  ylab("")+  
  theme_classic()+
  theme(axis.ticks.y = element_blank(), 
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.x = element_line(linetype = 1, size = 0.5),
        axis.ticks.x = element_line(linetype = 1, size = 0.5),
        plot.title = element_text(size = 11, family = "sans", hjust = 0.5, face="bold"),
        axis.text = element_text(size = 10, family = "sans"))
Inf_Fig3

# Figure 5
ggarrange(Inf_Fig3, SSTPlot, BWOxPlot, nrow = 1)
