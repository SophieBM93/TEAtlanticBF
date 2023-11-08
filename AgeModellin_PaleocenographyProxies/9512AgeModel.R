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

# 2) Set working Directory 
setwd("C:/Users/Sofia Barragan Monti/Documents/PhDBFMultivariate/TEAtlanticBF/AgeModellin_PaleocenographyProxies") 

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
Age<-read.csv(file="C:/Users/Sofia Barragan Monti/Documents/PhDBFMultivariate/TEAtlanticBF/AgeModellin_PaleocenographyProxies/Bacon_runs/GeoB9512-5/GeoB9512-5.csv", 
                  header=TRUE, stringsAsFactors=FALSE,sep=",",na.strings="-", check.names = F)

#d18O of Cibicidoides from https://doi.org/10.1594/PANGAEA.962968 
d18OCib<-read.csv(file="GeoB9512-5_Iso.csv", header=TRUE, stringsAsFactors=FALSE,
                  sep=",",na.strings="-", check.names = F)
d18OCib<-d18OCib[1:106,]

#d18O of Uvigerina from https://doi.org/10.1594/PANGAEA.962968 
d18OUvig<-read.csv(file="GeoB9512-5_iso.csv", header=TRUE, stringsAsFactors=FALSE,
                   sep=",",na.strings="-", check.names = F)
d18OUvig<-d18OUvig[107:133,]

#Bottom water temperatures
## Uvigerina Roberts et al. 2016 from (pending, not submitted yet)
BWTUvi<-read.csv(file="GeoB9512_BWT_Uvig(Rob.).csv", 
                 header=TRUE, stringsAsFactors=FALSE,sep=",",na.strings="-", check.names = F)
colnames(BWTUvi)<- c("Depth [cm]", "BWT[deg-C]")

#Ice Effect from doi.org/10.1038/35089060 
iceEffect<-read.csv(file="GlobalD18O.csv", 
                    header=TRUE, stringsAsFactors=FALSE,sep=",",na.strings="-", check.names = F)

########################################################
# 6) Run age model & save age model data
#Reload old run
Bacon("GeoB9512-5",ask=FALSE,ssize=1.5*nmc,d.min=d.min,d.max=d.max,acc.mean=50, 
      mem.mean=0.7,mem.strength=5,acc.shape=1.5, run=FALSE)
# Run Bacon
Bacon(core="GeoB9512-5", ask=FALSE,ssize=1.5*nmc,d.min=d.min,
      d.max=d.max,acc.mean=50,d.by = 1, mem.mean=0.7,mem.strength=5,
      acc.shape=1.5,t.a=20,t.b=21, res=20)

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
write.csv(dfAllChron,file=sprintf("OPGeoB9512-5_AgeModel.csv"), row.names=FALSE)

########################################################
# 7) Write out all simulated ages for the proxy depths
#Isotopes Cibicidoides
d18OMCAG<-sapply(d18OCib[,1], Bacon.Age.d) # Extract ages with Bacon for Isotopes
d18Oneff<-dim(d18OMCAG)[1] # Determine number of columns
d18Ochron.dat<-t(d18OMCAG)[,(d18Oneff-nmc+1):d18Oneff] # Tailor to desired number of iterations (nmc)
d18Ochron.median<-apply(d18Ochron.dat,1,median) # Calculate median age

#Isotopes Uvigerina
d18OMCAGUvig<-sapply(d18OUvig[,1], Bacon.Age.d) # Extract ages with Bacon for Isotopes
d18OneffUvig<-dim(d18OMCAGUvig)[1] # Determine number of columns
d18Ochron.datUvig<-t(d18OMCAGUvig)[,(d18OneffUvig-nmc+1):d18OneffUvig] # Tailor to desired number of iterations (nmc)
d18Ochron.medianUvig<-apply(d18Ochron.datUvig,1,median) # Calculate median age

#Bottom Water Temperatures (BWT) Uvigerina
BWT_MCAG<-sapply(BWTUvi[,1], Bacon.Age.d) # Extract ages with Bacon for BWT_Uvi
BWTneff<-dim(BWT_MCAG)[1] # Determine number of columns
BWTchron.dat<-t(BWT_MCAG)[,(BWTneff-nmc+1):BWTneff] # Tailor to desired number of iterations (nmc)
BWTchron.median<-apply(BWTchron.dat,1,median) # Calculate median age

########################################################
# 8) Add uncertainty to proxies
#d18O Cibicidoides
d18Odat.unc<-matrix(0,nrow=length(d18OCib[,2]),ncol=nmc,byrow=FALSE) # Prepare matrix for d18O + noise
d18Odat.prx<-matrix(d18OCib[,2],nrow=length(d18OCib[,2]),ncol=nmc,byrow=FALSE) # Prepare matrix with data
d18Osig.prx<-matrix(0,nrow=length(d18OCib[,2]),ncol=nmc,byrow=FALSE) # Prepare matrix for noise
for(j in 1:length(d18OCib[,2])){d18Osig.prx[j,]<-rnorm(n=nmc,sd=0.06)} # Put noise into matrix

d18Odat.unc<-d18Odat.prx+d18Osig.prx # Add noise to data and store in matrix

#d13C Cibicidoides
d13Cdat.unc<-matrix(0,nrow=length(d18OCib[,3]),ncol=nmc,byrow=FALSE) # Prepare matrix for d13C + noise
d13Cdat.prx<-matrix(d18OCib[,3],nrow=length(d18OCib[,3]),ncol=nmc,byrow=FALSE) # Prepare matrix with data
d13Csig.prx<-matrix(0,nrow=length(d18OCib[,3]),ncol=nmc,byrow=FALSE) # Prepare matrix for noise
for(j in 1:length(d18OCib[,3])){d13Csig.prx[j,]<-rnorm(n=nmc,sd=0.03)} # Put noise into matrix

d13Cdat.unc<-d13Cdat.prx+d13Csig.prx # Add noise to data and store in matrix

#d18O Uvigerina
d18Odat.uncUvig<-matrix(0,nrow=length(d18OUvig[,2]),ncol=nmc,byrow=FALSE) # Prepare matrix for d18O + noise
d18Odat.prxUvig<-matrix(d18OUvig[,2],nrow=length(d18OUvig[,2]),ncol=nmc,byrow=FALSE) # Prepare matrix with data
d18Osig.prxUvig<-matrix(0,nrow=length(d18OUvig[,2]),ncol=nmc,byrow=FALSE) # Prepare matrix for noise
for(j in 1:length(d18OUvig[,2])){d18Osig.prxUvig[j,]<-rnorm(n=nmc,sd=0.06)} # Put noise into matrix

d18Odat.uncUvig<-d18Odat.prxUvig+d18Osig.prxUvig # Add noise to data and store in matrix

#d13C Uvigerina
d13Cdat.uncUvig<-matrix(0,nrow=length(d18OUvig[,3]),ncol=nmc,byrow=FALSE) # Prepare matrix for d13C + noise
d13Cdat.prxUvig<-matrix(d18OUvig[,3],nrow=length(d18OUvig[,3]),ncol=nmc,byrow=FALSE) # Prepare matrix with data
d13Csig.prxUvig<-matrix(0,nrow=length(d18OUvig[,3]),ncol=nmc,byrow=FALSE) # Prepare matrix for noise
for(j in 1:length(d18OUvig[,3])){d13Csig.prxUvig[j,]<-rnorm(n=nmc,sd=0.03)} # Put noise into matrix

d13Cdat.uncUvig<-d13Cdat.prxUvig+d13Csig.prxUvig # Add noise to data and store in matrix

#BWT from Uvigerina Mg/Ca
BWTdat.unc<-matrix(0,nrow=length(BWTUvi[,2]),ncol=nmc,byrow=FALSE) # Prepare matrix for BWT_Uvi + noise
BWTdat.prx<-matrix(BWTUvi[,2],nrow=length(BWTUvi[,2]),ncol=nmc,byrow=FALSE) # Prepare matrix with data
BWTsig.prx<-matrix(0,nrow=length(BWTUvi[,2]),ncol=nmc,byrow=FALSE) # Prepare matrix for noise
for(j in 1:length(BWTUvi[,2])){BWTsig.prx[j,]<-rnorm(n=nmc,sd=0.46)} # Put noise into matrix

BWTdat.unc<-BWTdat.prx+BWTsig.prx # Add noise to data and store in matrix

#Global Ice effect
IEdat.unc<-matrix(0,nrow=length(iceEffect[,2]),ncol=nmc,byrow=FALSE) # Prepare matrix for ice effect + noise
IEdat.prx<-matrix(iceEffect[,2],nrow=length(iceEffect[,2]),ncol=nmc,byrow=FALSE) # Prepare matrix with mean data
IEsig.prx<-matrix(0,nrow=length(iceEffect[,2]),ncol=nmc,byrow=FALSE) # Prepare matrix for noise
for(j in 1:length(iceEffect[,2])){IEsig.prx[j,]<-sample((iceEffect[j,3]*100):(iceEffect[j,4]*-100),nmc, replace=TRUE)/100} # Put noise into matrix

IEdat.unc<-IEdat.prx+IEsig.prx # Add noise to data and store in matrix

########################################################
# 9) Interpolate ensembles to the median age

#d18O of Cibicidoides 
d18OMedAge<-matrix(NA, nrow=length(BWTchron.median), ncol=nmc) # Matrix for interpolated values
for (i in 1:nmc)
{
d18OMedAge[,i]<-approx(d18Ochron.dat[,i],d18Odat.unc[,i],BWTchron.median)$y 
}
d18OMedAge.mean<-matrix(apply(d18OMedAge,1,mean,na.rm=TRUE)) # Calculate mean d18O

#BWT data
BWTMedAge<-matrix(NA, nrow=length(BWTchron.median), ncol=nmc) # Matrix for interpolated values BWT_Uvi
for (i in 1:nmc)
{
BWTMedAge[,i]<-approx(BWTchron.dat[,i],BWTdat.unc[,i],BWTchron.median)$y 
}
BWTMedAge.mean<-matrix(apply(BWTMedAge,1,mean,na.rm=TRUE)) # Calculate mean d18O

##Ice effect
IEMedAged18O<-matrix(NA, nrow=length(BWTchron.median), ncol=nmc) # Matrix for interpolated values
for (i in 1:nmc)
{
IEMedAged18O[,i]<-approx(iceEffect[,1]*1000,IEdat.unc[,i],BWTchron.median)$y 
}

########################################################
# 10) Interpolate d13C/d18O Ensemble to the median age

#d18O - Cib
d18OMedAgeHigh<-matrix(NA, nrow=length(d18Ochron.median), ncol=nmc) # Matrix for interpolated values
for (i in 1:nmc)
{
d18OMedAgeHigh[,i]<-approx(d18Ochron.dat[,i],d18Odat.unc[,i],d18Ochron.median)$y 
}

d18OMedAgeHigh.mean<-matrix(apply(d18OMedAgeHigh,1,mean,na.rm=TRUE)) # Calculate mean d18O

#d13C - Cib
d13CMedAgeHigh<-matrix(NA, nrow=length(d18Ochron.median), ncol=nmc) # Matrix for interpolated values
for (i in 1:nmc)
{
d13CMedAgeHigh[,i]<-approx(d18Ochron.dat[,i],d13Cdat.unc[,i],d18Ochron.median)$y #we use d18O time scale (same depth)
}

d13CMedAgeHigh.mean<-matrix(apply(d13CMedAgeHigh,1,mean,na.rm=TRUE)) # Calculate mean d18O

#d18O Uvig
d18OMedAgeHighUvig<-matrix(NA, nrow=length(d18Ochron.medianUvig), ncol=nmc) # Matrix for interpolated values
for (i in 1:nmc)
{
  d18OMedAgeHighUvig[,i]<-approx(d18Ochron.datUvig[,i],d18Odat.uncUvig[,i],d18Ochron.medianUvig)$y 
}

d18OMedAgeHigh.meanUvig<-matrix(apply(d18OMedAgeHighUvig,1,mean,na.rm=TRUE)) # Calculate mean d18O

#d13C Uvig
d13CMedAgeHighUvig<-matrix(NA, nrow=length(d18Ochron.medianUvig), ncol=nmc) # Matrix for interpolated values
for (i in 1:nmc)
{
  d13CMedAgeHighUvig[,i]<-approx(d18Ochron.datUvig[,i],d13Cdat.uncUvig[,i],d18Ochron.medianUvig)$y #we use d18O time scale (same depth)
}

d13CMedAgeHigh.meanUvig<-matrix(apply(d13CMedAgeHighUvig,1,mean,na.rm=TRUE)) # Calculate mean d18O

########################################################
# 11) Calculate seawater d18O for ensemble

d18Ow <- ((BWTMedAge-16.9+(4*d18OMedAge))/4)+0.64             #d18Osw using Schackelton, 1984 +0.64 correction added
d18Ow.mean<-matrix(apply(d18Ow,1,mean,na.rm=TRUE))

########################################################
# 12) Extract quantiles from time series
probs <- c(0.025,0.16,0.84,0.975)
Qts18O <- t(apply(d18OMedAgeHigh,1,quantile,probs=probs,na.rm=TRUE))
Qts13C <- t(apply(d13CMedAgeHigh,1,quantile,probs=probs,na.rm=TRUE))
Qts18OUvig <- t(apply(d18OMedAgeHighUvig,1,quantile,probs=probs,na.rm=TRUE))
Qts13CUvig <- t(apply(d13CMedAgeHighUvig,1,quantile,probs=probs,na.rm=TRUE))
QtsBWT <- t(apply(BWTMedAge,1,quantile,probs=probs,na.rm=TRUE))
Qts18Ow <- t(apply(d18Ow,1,quantile,probs=probs,na.rm=TRUE))

########################################################
# 13) Save data 
dfd18O<-data.frame(cbind(d18OCib[,1],d18Ochron.median/1000,d18OMedAgeHigh.mean,Qts18O)) # Merge fields
colnames(dfd18O)<-c("Depth [cm]","Median Age [kyrs]","Mean d18O [o/oo]","Q0.025","Q0.16", "Q0.84", "Q0.975")
write.csv(dfd18O,file=sprintf("OPGeoB9512-5_Cib18OQts.csv"), row.names=FALSE)

dfd13C<-data.frame(cbind(d18OCib[,1],d18Ochron.median/1000,d13CMedAgeHigh.mean,Qts13C)) # Merge fields
colnames(dfd13C)<-c("Depth [cm]","Median Age [kyrs]","Mean d13C [o/oo]","Q0.025","Q0.16", "Q0.84", "Q0.975")
write.csv(dfd13C,file=sprintf("OPGeoB9512-5_Cib13CQts.csv"), row.names=FALSE)

dfd18OUvig<-data.frame(cbind(d18OUvig[,1],d18Ochron.medianUvig/1000,d18OMedAgeHigh.meanUvig,Qts18OUvig)) # Merge fields
colnames(dfd18OUvig)<-c("Depth [cm]","Median Age [kyrs]","Mean d18O [o/oo]","Q0.025","Q0.16", "Q0.84", "Q0.975")
write.csv(dfd18OUvig,file=sprintf("OPGeoB9512-5_Uvig18OQts.csv"), row.names=FALSE)

dfd13CUvig<-data.frame(cbind(d18OUvig[,1],d18Ochron.medianUvig/1000,d13CMedAgeHigh.meanUvig,Qts13CUvig)) # Merge fields
colnames(dfd13CUvig)<-c("Depth [cm]","Median Age [kyrs]","Mean d13C [o/oo]","Q0.025","Q0.16", "Q0.84", "Q0.975")
write.csv(dfd13CUvig,file=sprintf("OPGeoB9512-5_Uvig13CQts.csv"), row.names=FALSE)

dfBWT<-data.frame(cbind(BWTUvi[,1],BWTchron.median/1000,BWTMedAge.mean,QtsBWT)) # Merge fields
colnames(dfBWT)<-c("Depth [cm]","Median Age [kyrs]","BWT [deg C]","Q0.025","Q0.16", "Q0.84", "Q0.975")
write.csv(dfBWT,file=sprintf("OPGeoB9512-5_BWTUvigerina.csv"), row.names=FALSE)

dfdw<-data.frame(cbind(BWTchron.median/1000,d18Ow.mean,Qts18Ow)) # Merge fields
colnames(dfdw)<-c("Median Age [kyrs]","dw [permil]","Q0.025","Q0.16", "Q0.84", "Q0.975")
write.csv(dfdw,file=sprintf("OPGeoB9512-5_dw.csv"), row.names=FALSE)

################################
# 14) XRF to calibrate Heinrich stadials and Younger Dryas https://doi.org/10.1029/2018PA003359
XRF<-read.csv(file="GeoB9512-5_XRF.csv", 
              header=TRUE, stringsAsFactors=FALSE,sep=",",na.strings="-", check.names = F)

# Write out all simulated ages for depths 

XRFMCAG<-sapply(XRF[,1], Bacon.Age.d) # Extract ages with Bacon
XRFneff<-dim(XRFMCAG)[1] # Determine number of columns
XRFchron.dat<-t(XRFMCAG)[,(XRFneff-nmc+1):XRFneff] # Tailor to desired number of iterations (nmc)
XRFchron.median<-apply(XRFchron.dat,1,median) # Calculate median age

XRF_Age<-data.frame(cbind(XRFchron.median/1000,XRF)) # Merge fields
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

dfd18O<-read.csv("OPGeoB9512-5_Cib18OQts.csv", check.names = F)
dfd13C<-read.csv("OPGeoB9512-5_Cib13CQts.csv", check.names=F)

dfd18O3<-read.csv("OPGeoB9512-5_Uvig18OQts.csv", check.names = F)
dfd13C3<-read.csv("OPGeoB9512-5_Uvig13CQts.csv", check.names=F)

dfBWT<-read.csv("OPGeoB9512-5_BWTUvigerina.csv", check.names=F)

dfdw<-read.csv("OPGeoB9512-5_dw.csv", check.names=FALSE)

XRF <- read.csv("OPGeoB9512-5_XRF_Age.csv", header = T, check.names = F)

Globald18O <- read.csv(file="GlobalD18O.csv", check.names = F, header=TRUE, 
                       stringsAsFactors=FALSE,sep=",",na.strings="-")

########################################################
# 2) Define Key Climatic Periods (As reported in the Manuscript)

YD <- tibble(ymin = 14.7, ymax =13.2, xmin = -Inf, xmax = Inf)
HS1 <- tibble(ymin = 18.2, ymax =15.8, xmin = -Inf, xmax = Inf)
HS2 <- tibble(ymin = 25.84, ymax =24.29, xmin = -Inf, xmax = Inf)
LGM <- tibble(ymin = 23, ymax =19, xmin = -Inf, xmax = Inf)

########################################################
# 3) Plot Elemental ratios (uncalibrated) with the new age model to confirm HS2, HS1 and YD
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
# 4) Plot Paleocenographic proxies with uncertainties

## d18O raw data
d18Oplotraw<- ggplot(dfd18O, aes(`Mean d18O [o/oo]`, `Median Age [kyrs]`))+
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
  geom_ribbon(dfd18O, mapping=aes(xmin = `Q0.16`, xmax = `Q0.84`, y=`Median Age [kyrs]`),
              inherit.aes = F, fill = "gray", alpha=0.3)+
  geom_lineh(dfd18O, mapping=aes(`Mean d18O [o/oo]`, `Median Age [kyrs]`),
             colour="black", linetype = 2, linewidth = 0.25)+
  geom_point(dfd18O, mapping=aes(`Mean d18O [o/oo]`, `Median Age [kyrs]`),
             colour="black", size = 1, alpha = 0.8)+
  scale_y_reverse(guide="prism_minor",
                  breaks = seq(50,0,by=-5),
                  limits = c(27,0),
                  expand = c(0,0),
                  minor_breaks=seq(50, 0, -1))+
  scale_x_reverse(guide="prism_minor",
                  breaks = seq(20,-20,by=-1),
                  limits = c(4.5, 1.5),
                  expand = c(0,0),
                  minor_breaks=seq(20,-20, -0.5))+
  labs(title="d18OCib  [%. VSMOW]")+
  xlab("")+
  ylab("Age (kyrs BP)")+
  theme_classic()+
  theme(axis.title.y = element_text(size = 11, family = "sans", hjust = 0.5, face="bold"),
        axis.text.y = element_text(size = 10, family = "sans"),
        axis.line.x = element_line(linetype = 1, size = 0.5),
        axis.ticks.x = element_line(linetype = 1, size = 0.5),
        plot.title = element_text(size = 11, family = "sans", hjust = 0.5, face="bold"),
        axis.text.x = element_text(size = 10, family = "sans"))
d18Oplotraw

## d18Osw
d18Oplotsw <- ggplot(dfdw, aes(`dw [permil]`, `Median Age [kyrs]`))+
  geom_vline(xintercept = 0.54, linetype = 2, linewidth = 0.5, alpha=0.5)+
  geom_ribbon(dfdw, mapping=aes(xmin = `Q0.16`, xmax = `Q0.84`, y=`Median Age [kyrs]`),
              inherit.aes = F, fill = "#7d9ec2", alpha=0.3)+ 
  geom_lineh(dfdw, mapping=aes(`dw [permil]`, `Median Age [kyrs]`), 
             colour="#7d9ec2", linetype = 2, linewidth = 0.25)+
  geom_point(dfdw, mapping=aes(`dw [permil]`, `Median Age [kyrs]`), 
             colour="blue", size=1, alpha=0.5)+
  geom_lineh(Globald18O, mapping= aes(`delta O-18`, `Age [ka]`), 
             colour="black", linetype = 1, linewidth = 0.25, alpha=0.6)+
  scale_y_reverse(guide="prism_minor",
                  breaks = seq(50,0,by=-5),
                  limits = c(27,0),
                  expand = c(0,0),
                  minor_breaks=seq(50, 0, -1))+
  scale_x_continuous(guide="prism_minor",
                     breaks = seq(-20,20,by=1),
                     limits = c(-1, 3),
                     expand = c(0,0),
                     minor_breaks=seq(-20, 20, 0.5))+
  labs(title="d18Osw  [%. VSMOW]")+
  xlab("")+
  ylab("Age (kyrs BP)")+
  theme_classic()+
  theme(axis.line.y =element_blank(), 
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.x = element_line(linetype = 1, size = 0.5),
        axis.ticks.x = element_line(linetype = 1, size = 0.5),
        plot.title = element_text(size = 11, family = "sans", hjust = 0.5, face="bold"),
        axis.text.x = element_text(size = 10, family = "sans"))
d18Oplotsw

## d13C
d13Cplot <- ggplot(dfd13C, aes(`Mean d13C [o/oo]` , `Median Age [kyrs]`))+
  geom_hline(yintercept = 6.6, linetype = 2, linewidth = 0.5, alpha=0.5)+
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            data = LGM, alpha = 0.6, fill = "#c7c7c7", inherit.aes = FALSE)+
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            data = HS1, alpha = 0.4, fill = "#A6B578",inherit.aes = FALSE)+
  geom_hline(yintercept = 16.59, linetype=2, linewidth=0.5, alpha=0.5)+
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            data = HS2, alpha = 0.4, fill = "#A6B578",inherit.aes = FALSE)+
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            data = YD, alpha = 0.6, fill = "#c7c7c7", inherit.aes = FALSE)+
  geom_vline(xintercept = 0.53, linetype = 2, size = 0.5, alpha=0.5)+
  geom_ribbon(dfd13C, mapping=aes(xmin = `Q0.16`, xmax = `Q0.84`), fill = "#A6B578", alpha=0.3)+
  geom_lineh(dfd13C, mapping=aes(`Mean d13C [o/oo]` , `Median Age [kyrs]`),
             colour="darkgreen", linetype = 2, size = 0.25)+
  geom_point(dfd13C, mapping=aes(`Mean d13C [o/oo]` , `Median Age [kyrs]`),
             colour="#A6B578", size = 1, alpha = 0.5)+
  scale_y_reverse(guide="prism_minor",
                  breaks = seq(50,0,by=-5),
                  limits = c(27,0),
                  expand = c(0,0),
                  minor_breaks=seq(50, 0, -1))+
  scale_x_continuous(guide="prism_minor",
                     breaks = seq(-2,20,by=0.5),
                     limits = c(-0.8,2),
                     expand = c(0,0),
                     minor_breaks=seq(-2,2, 0.25))+
  labs(title="d13C [%. VPDB]")+ xlab("")+ ylab("Age (Kyrs)")+
  theme_classic()+
  theme(axis.line.y =element_blank(), 
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.x = element_line(linetype = 1, size = 0.5),
        axis.ticks.x = element_line(linetype = 1, size = 0.5),
        plot.title = element_text(size = 11, family = "sans", hjust = 0.5, face="bold"),
        axis.text.x = element_text(size = 10, family = "sans"))
d13Cplot

## Plot all proxies together
plot_grid(d18Oplotraw, d18Oplotsw, d13Cplot, align="hv",
          nrow = 1, hjust = 0, vjust = 0)

## 5) Supplementary Figure 1
# Import Environmental Data
ea1 <- read.csv("9512_EA.csv", check.names = F)
ea1 <- replace(ea1,is.na(ea1),0)                    #Replacing NA values by 0
ea1 <- ea1[rowSums(ea1[,26:27])>0,]                 #Delete 0 rows



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
TestPlot                                         #Black are porcelaneous, dark gray Calcareous

#Stress Species
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
  theme_classic()+
  theme(axis.ticks.y = element_blank(), 
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.x = element_line(linetype = 1, size = 0.5),
        axis.ticks.x = element_line(linetype = 1, size = 0.5),
        plot.title = element_text(size = 11, family = "sans", hjust = 0.5, face="bold"),
        axis.text = element_text(size = 10, family = "sans"))
SSPlot

#Diveristy plots
dPlot <- ggplot(ea1, aes(`Fisher Alpha Index`, `Age [ka]`))+
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
  theme(axis.ticks.y = element_blank(), 
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.x = element_line(linetype = 1, size = 0.5),
        axis.ticks.x = element_line(linetype = 1, size = 0.5),
        plot.title = element_text(size = 11, family = "sans", hjust = 0.5, face="bold"),
        axis.text = element_text(size = 10, family = "sans"))
dPlot

#Most Abundant species
# install.packages("reshape")
library("reshape")

ea1_b <- cbind(ea1$`Age [ka]`, ea1[,26:42])

melted <- melt(ea1_b, id.vars = "ea1$`Age [ka]`")
str(melted)
names(melted)[1] <- "Age"


AbundSp <- ggplot(melted, aes(x=value,y= Age, fill=variable))+
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
  scale_fill_manual(values = c("grey","grey","grey","grey","black","black","black",
                    "black", "black","black", "black", "black", "black", "black",
                    "black", "grey","grey","grey","grey"))+
  theme_classic()+
  theme(legend.position = "none")+
  facet_wrap(~ variable, nrow=1)

AbundSp 

#Supplementary Figure 1
ggarrange(AbundSp ,d18Oplotraw, d18Oplotsw, TestPlot, SSPlot, dPlot, 
          nrow = 1, widths = c(17,3,3,3,3,2.5))


## 6) Figure 3
## Data from KNR166-2-29JPC extracted from https://doi.org/10.1594/PANGAEA.936747
# age model was calculated with the same methodology 
# using radiocarbon ages reported https://doi.org/10.1594/PANGAEA.936747 (https://doi.pangaea.de/10.1594/PANGAEA.831673)

NW_Atl <- read.csv("KNR166-2-29JPC_Cib13CQts.csv", check.names = F)

CibSpPlot <- ggplot(ea1, aes(`Cibicidoides spp.`, `Age [ka]`))+
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
  geom_areah(fill = "grey70")+
  geom_vline(xintercept = 15, linetype = 2, linewidth = 0.5, alpha=0.5)+
  scale_y_reverse(guide="prism_minor",
                  breaks = seq(50,0,by=-5),
                  limits = c(27,0),
                  expand = c(0,0),
                  minor_breaks=seq(50, 0, -1))+
  scale_x_continuous(guide="prism_minor",
                     breaks = seq(-100,200,by=20),
                     limits = c(0, 50),
                     expand = c(0,0),
                     minor_breaks=seq(-100,200,by=10))+
  labs(title="Cibicidoides spp")+ xlab("")+ ylab("Age (Kyrs)")+
  theme_classic()+
  theme(axis.line.x = element_line(linetype = 1, size = 0.5),
        axis.ticks.x = element_line(linetype = 1, size = 0.5),
        plot.title = element_text(size = 11, family = "sans", hjust = 0.5, face="bold"),
        axis.text.x = element_text(size = 10, family = "sans"))
CibSpPlot

d13Cplot2 <- ggplot(dfd13C, aes(`Mean d13C [o/oo]` , `Median Age [kyrs]`))+
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
  geom_vline(xintercept = 0.53, linetype = 2, size = 0.5, alpha=0.5)+
  geom_ribbon(dfd13C, mapping=aes(xmin = `Q0.16`, xmax = `Q0.84`), fill = "#A6B578", alpha=0.3)+
  geom_lineh(dfd13C, mapping=aes(`Mean d13C [o/oo]` , `Median Age [kyrs]`), 
             colour="#A6B578", linetype = 2, size = 0.25)+
  geom_point(dfd13C, mapping=aes(`Mean d13C [o/oo]` , `Median Age [kyrs]`), 
             colour="#A6B578", size = 1, alpha = 0.5)+
  geom_ribbon(NW_Atl, mapping=aes(xmin = `Q0.16`, xmax = `Q0.84`), fill = "#993399", alpha=0.3)+
  geom_lineh(NW_Atl, mapping=aes(`Mean d13C [o/oo]` , `Median Age [kyrs]`), 
             colour="#993399", linetype = 2, size = 0.25)+
  geom_point(NW_Atl, mapping=aes(`Mean d13C [o/oo]` , `Median Age [kyrs]`), 
             colour="#993399", size = 1, alpha = 0.5)+
  geom_lineh(ea1, mapping=aes(`Dd13C (9512-29JPC)` , `Average Age [kyrs] Dd13C`),
             colour="black", linetype = 2, size = 0.25)+
  geom_point(ea1, mapping=aes(`Dd13C (9512-29JPC)` , `Average Age [kyrs] Dd13C`),
             colour="black", size = 1, alpha = 0.5)+
  scale_y_reverse(guide="prism_minor",
                  breaks = seq(50,0,by=-5),
                  limits = c(27,0),
                  expand = c(0,0),
                  minor_breaks=seq(50, 0, -1))+
  scale_x_continuous(guide="prism_minor",
                     breaks = seq(-2,20,by=0.5),
                     limits = c(-2,2),
                     expand = c(0,0),
                     minor_breaks=seq(-2,2, 0.25))+
  labs(title="d13C [%. VPDB]")+ xlab("")+ ylab("")+
  theme_classic()+
  theme(axis.ticks.y = element_blank(), 
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.x = element_line(linetype = 1, size = 0.5),
        axis.ticks.x = element_line(linetype = 1, size = 0.5),
        plot.title = element_text(size = 11, family = "sans", hjust = 0.5, face="bold"),
        axis.text.x = element_text(size = 10, family = "sans"))
d13Cplot2

## Dissolve Oxygen concentration from BFOI
# 2) Oxic zones (https://www.nature.com/articles/s41598-022-05295-8)

HighOx <- tibble(xmin = 3, xmax =6, ymin = -Inf, ymax = Inf)
LowOx <- tibble(xmin = 1.5, xmax =3, ymin = -Inf, ymax = Inf)
Subox <- tibble(xmin = 0.3, xmax =1.5, ymin = -Inf, ymax = Inf)
Dysox <- tibble(xmin = 0, xmax =0.3, ymin = -Inf, ymax = Inf)

BWOxPlot <- ggplot(ea1, aes(`Inferred Dissolve Oxygen Concentration [mL/L]`, `Age [ka]`))+      #Plots the enhanced BFOI
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
            data = LowOx, alpha = 0.6, fill = "#c7c7c7", inherit.aes = FALSE)+
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
            data = Dysox, alpha = 0.6, fill = "#c7c7c7",inherit.aes = FALSE)+
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
  theme_classic()+
  theme(axis.ticks.y = element_blank(), 
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.x = element_line(linetype = 1, size = 0.5),
        axis.ticks.x = element_line(linetype = 1, size = 0.5),
        plot.title = element_text(size = 11, family = "sans", hjust = 0.5, face="bold"),
        axis.text.x = element_text(size = 10, family = "sans"))
BWOxPlot

#TROX Model
ModOM <- tibble(xmin = 33, xmax =66, ymin = -Inf, ymax = Inf)

TROXPlot <- ggplot(ea1, aes(`Infaunals[%]`, `Age [ka]`)) +
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax), 
            data = ModOM, alpha = 0.6, fill = "#c7c7c7", inherit.aes = FALSE)+
  geom_lineh(linewidth=0.5, colour="black", linetype=1) +
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
  labs(title="Infaunal (%)") +
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
TROXPlot

#XRF Mn
MnPlot <- ggplot(XRF_Age, aes(Mn, `Median Age [kyrs]`)) +
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = 30, xmax = 50), 
            data = ModOM, alpha = 0.6, fill = "#c7c7c7", inherit.aes = FALSE)+
  geom_lineh(linewidth=0.5, colour="black", linetype=1) +
  scale_y_reverse(guide="prism_minor",
                  breaks = seq(50,0,by=-5),
                  limits = c(27,0),
                  expand = c(0,0),
                  minor_breaks=seq(50, 0, -1))+
  scale_x_continuous(guide="prism_minor",
                     breaks = seq(0,100,by=20),
                     limits = c(0,60),
                     expand = c(0,0),
                     minor_breaks=seq(0,100, 10))+
  labs(title="XRF-Mn") +
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

MnPlot

#XRF TiCa
TiCaPlot <- ggplot(XRF_Age, aes(`Ti/Ca`, `Median Age [kyrs]`)) +
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = 0.1, xmax = 0.6), 
            data = ModOM, alpha = 0.6, fill = "#c7c7c7", inherit.aes = FALSE)+
  geom_lineh(linewidth=0.5, colour="black", linetype=1) +
  scale_y_reverse(guide="prism_minor",
                  breaks = seq(50,0,by=-5),
                  limits = c(27,0),
                  expand = c(0,0),
                  minor_breaks=seq(50, 0, -1))+
  scale_x_continuous(guide="prism_minor",
                     breaks = seq(0,100,by=0.2),
                     limits = c(0,0.6),
                     expand = c(0,0),
                     minor_breaks=seq(0,100, 0.1))+
  labs(title="XRF-Mn") +
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
TiCaPlot

#XRF FeCa
FeCaPlot <- ggplot(XRF_Age, aes(`Fe/Ca`, `Median Age [kyrs]`)) +
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = 2, xmax = 8.5), 
            data = ModOM, alpha = 0.6, fill = "#c7c7c7", inherit.aes = FALSE)+
  geom_lineh(linewidth=0.5, colour="black", linetype=1) +
  scale_y_reverse(guide="prism_minor",
                  breaks = seq(50,0,by=-5),
                  limits = c(27,0),
                  expand = c(0,0),
                  minor_breaks=seq(50, 0, -1))+
  scale_x_continuous(guide="prism_minor",
                     breaks = seq(0,100,by=2),
                     limits = c(0,10),
                     expand = c(0,0),
                     minor_breaks=seq(0,100, 1))+
  labs(title="XRF-Mn") +
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

FeCaPlot

#Figure 3
ggarrange(CibSpPlot, d13Cplot2, BWOxPlot, TROXPlot, 
          dPlot, MnPlot, TiCaPlot, FeCaPlot,
          nrow = 1, widths = c(4,6,4,4,4,3,3,3))


## 7) Figure 6
BFDRPlot <- ggplot(ea1, aes(BFDR, `Average Age [kyrs] BFDR`))+      
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = -10, xmax = 0), 
            data = LowOx, alpha = 0.6, fill = "#c7c7c7", inherit.aes = FALSE)+
  scale_y_reverse(guide="prism_minor",
                  breaks = seq(50,0,by=-5),
                  limits = c(27,0),
                  expand = c(0,0),
                  minor_breaks=seq(50, 0, -1))+
  scale_x_continuous(guide="prism_minor",
                     breaks = seq(-100,100,by=10),
                     limits = c(-30,30),
                     expand = c(0,0),
                     minor_breaks=seq(-100,100, 5))+
  geom_lineh()+
  labs(title="Benthic Foraminifera Diversification Rates") +
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
BFDRPlot

dPlot2 <- ggplot(ea1, aes(`Fisher Alpha Index`, `Age [ka]`))+
  geom_rect(mapping = aes(ymin = ymin, ymax = ymax, xmin = 7.9, xmax = 12.99), 
            data = LowOx, alpha = 0.6, fill = "#c7c7c7", inherit.aes = FALSE)+
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
  theme(axis.line = element_line(linetype = 1, size = 0.5),
        axis.ticks = element_line(linetype = 1, size = 0.5),
        plot.title = element_text(size = 11, family = "sans", hjust = 0.5, face="bold"),
        axis.text = element_text(size = 10, family = "sans"))
dPlot2


#Figure 3
ggarrange(dPlot2, BFDRPlot, BWOxPlot, TROXPlot,
          nrow = 1, widths = c(3,4,3,3))


#####################################
# Part 3-Optional: Plots Global Variables
# 1) Upload data (sources cited) 

##Global CO2 https://doi.pangaea.de/10.1594/PANGAEA.871265
GlobalCO2 <- read.csv(file="Kohler et al., 2017_CO2.csv", check.names = F, 
                      header=TRUE, stringsAsFactors=FALSE,sep=",",na.strings="-")

##Global Mean Ocean Temperature https://www.nature.com/articles/s41586-021-03984-4
dGMST <- read.csv("gmst_Osmanetal.,2021.csv", check.names = F)
dGMST$`Age min (kyrs BP)`<-dGMST$`Age min (kyrs BP)`/1000

##Ph/Th from the Atlantic Ocean https://www.nature.com/articles/nature02494
PhTh_A <- read.csv(file="PhTh_McManus.csv", check.names = F, header=TRUE, 
                   stringsAsFactors=FALSE,sep=",",na.strings="-")
colnames(PhTh_A)<-c("Depth [cm]", "Age [kyrs]", "Pa/Th238", "Error238", "Pa/Th232","Error232")

#https://www.nature.com/articles/nature14059
PhTh_C <- read.csv(file="Bohm2015.csv", check.names = F, header=TRUE, 
                   stringsAsFactors=FALSE,sep=",",na.strings="-")
colnames(PhTh_C)<-c("Event", "Depth [cm]", "Age [kyrs]", "231Pa/230Th", "2sd %", "2SD")

########################
# 2) Plots the variables

#Atlantic Ph/Th
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
  geom_lineh(PhTh_A, mapping=aes(`Pa/Th238`, `Age [kyrs]`),
             colour="#9999CC", linetype = 1, size = 0.25)+
  geom_errorbar(PhTh_A,mapping=aes(y= `Age [kyrs]`,xmin = `Pa/Th238` - Error238,
                                   xmax = `Pa/Th238` + Error238),
                colour= "#9999CC", linetype = 1, width = 0.5)+
  geom_lineh(PhTh_A, mapping=aes(`Pa/Th232`, `Age [kyrs]`),
             colour="#66CC99", linetype = 1, size = 0.25)+
  geom_errorbar(PhTh_A,mapping=aes(y= `Age [kyrs]`,xmin = `Pa/Th232` - Error232,
                                   xmax = `Pa/Th232` + Error232),
                colour= "#66CC99", linetype = 1, width = 0.5)+
  geom_lineh(PhTh_C, mapping=aes(`231Pa/230Th`, `Age [kyrs]`),
             colour="#CC6666", linetype = 1, size = 0.25)+
  geom_errorbar(PhTh_C,mapping=aes(y= `Age [kyrs]`,xmin = `231Pa/230Th` - `2SD`,
                                   xmax = `231Pa/230Th` + `2SD`),
                colour="#CC6666", linetype = 1, width = 0.1, size=0.1)+
  scale_y_reverse(guide="prism_minor",
                  breaks = seq(180,0,by=-5),
                  limits = c(27,0),
                  expand = c(0,0),
                  minor_breaks=seq(180, 0, -1))+
  scale_x_reverse(guide="prism_minor",
                  breaks = seq(0.14,0,by=-0.02),
                  limits = c(0.14,0),
                  expand = c(0,0),
                  minor_breaks=seq(0.14,0,-0.01))+
  labs(title="231Pa/230Th")+
  xlab("")+ ylab("")+ theme_classic()+
  theme(axis.line.y =element_blank(), 
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
    plot.title = element_text(size = 14, family = "sans", hjust = 0.5),
        axis.text.x = element_text(size = 10, family = "sans"), 
        axis.line.x = element_line(linetype = 1, size = 0.5),
        axis.ticks.x = element_line(linetype = 1, size = 0.5))
PhThPlot

# pdf(file = "GeoB9512-5_PaleocForDiv_Nov23.pdf", paper = "a4", bg="transparent", width = 10, height = 5)
ggarrange(d18Oplotraw, d18Oplotsw, d13Cplot, PhThPlot, nrow = 1)
# dev.off()

#Global CO2
GlobalCO2Plot <- ggplot()+
  geom_lineh(GlobalCO2, mapping= aes(`CO2 [mmol/mol]`, `Age [ka BP]`), 
             colour="#B08267", linetype = 1, size = 0.4)+
  geom_point(GlobalCO2, mapping= aes(`CO2 std dev [±]`, `Age [ka BP]`), 
             colour="#B08267", size = 0.4)+
  scale_y_reverse(guide="prism_minor",
                  breaks = seq(50,0,by=-5),
                  limits = c(27,0),
                  expand = c(0,0),
                  minor_breaks=seq(50, 0, -1))+
  scale_x_continuous(guide="prism_minor",
                     breaks = seq(150,350,by=50),
                     limits = c(150,350),
                     expand = c(0,0),
                     minor_breaks=seq(150,350, 25))+
  labs(title="Global CO2 [mmol/mol]")+
  xlab("")+ 
  ylab("Age (Kyrs)")+
  theme_classic()+
  theme(axis.line.y =element_blank(), 
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
    plot.title = element_text(size = 14, family = "sans", hjust = 0.5),
        axis.text.x = element_text(size = 10, family = "sans"), 
        axis.line.x = element_line(linetype = 1, size = 0.5),
        axis.ticks.x = element_line(linetype = 1, size = 0.5))
GlobalCO2Plot


#Define colors and plot Global Mean Seasurface Temprature (GMST)
my_colors <- RColorBrewer::brewer.pal(6, "Blues_Pastel")[2:6]
gmst1 <- ggplot(dGMST, aes(`50th`, `Age min (kyrs BP)`))+
  geom_ribbon(dGMST,mapping=aes(xmin=`5th`, xmax=`95th`), fill = "#6A75BA", alpha=0.7)+
  geom_ribbon(dGMST,mapping=aes(xmin=`10th`, xmax=`90th`), fill = "#7A86CC", alpha=0.7)+
  geom_ribbon(dGMST,mapping=aes(xmin=`20th`, xmax=`80th`), fill = "#929CDE", alpha=0.7)+
  geom_ribbon(dGMST,mapping=aes(xmin=`30th`, xmax=`70th`), fill = "#B4BBFA", alpha=0.7)+
  geom_ribbon(dGMST,mapping=aes(xmin=`40th`, xmax=`60th`), fill = "#D3EEFF", alpha=0.7)+
  geom_lineh(linetype=2, size=0.5, colour="#7B60BD")+
  scale_y_reverse(guide="prism_minor",
                  breaks = seq(50,0,by=-5),
                  limits = c(27,0),
                  expand = c(0,0),
                  minor_breaks=seq(50, 0, -1))+
  scale_x_continuous(guide="prism_minor",
                     breaks = seq(-10,2,by=2),
                     limits = c(-10,2),
                     expand = c(0,0),
                     minor_breaks=seq(-10, 2, 1))+
  labs(title="dGMST [?C]")+
  xlab("")+ 
  ylab("Age (kyrs)")+
  theme_classic()+
  theme(axis.line.y =element_blank(), 
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
       plot.title = element_text(size = 11, family = "sans", hjust = 0.5, face="bold"),
        axis.text.x = element_text(size = 10, family = "sans"))
gmst1

GlobalVarib <- plot_grid(d18Oplotsw2,dfBWTplot,GlobalCO2Plot, dfd13Cplot, PhThPlot, gmst1, align="hv",
                          nrow = 1, hjust = 0, vjust = 0)
GlobalVarib

pdf(file = "Global.pdf", paper = "a4r", bg="transparent", width = 10, height = 5)
GlobalVarib
dev.off()
