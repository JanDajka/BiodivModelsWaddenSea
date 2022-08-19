setwd(".../Data/Macrozoobenthos")

######################################################################

library(lattice)
library(dplyr)
library(mgcv)
library(piecewiseSEM)
library(lme4)
library(psych)
library(plyr)
library(tidyverse)
library(lubridate)
library(vegan)
library(sjPlot)
library(MASS)
library(nlme)
library(sjPlot)
library(reshape2)
library(gdm)

###################################################################### Standing Diversity SEM

mzbwadden <- read.csv("MARISCO_mzb_wadden.csv")

mzbwadden$StationID <- as.factor(mzbwadden$StationID)

#Are the missing values?

levels(mzbwadden$StationID)
colSums(is.na(mzbwadden))
dim(mzbwadden)

mzbwadden <- na.exclude(mzbwadden) 

dim(mzbwadden)

levels(mzbwadden$StationID)
summary(mzbwadden)
# Outliers

par(mfrow = c(1, 3))
boxplot(mzbwadden$Sy, 
        main = "Species richness y")
hist(mzbwadden$Sy)
dotchart(mzbwadden$Sy)
boxplot(mzbwadden$S, 
        main = "Species richness a")
hist(mzbwadden$S)
dotchart(mzbwadden$S)
boxplot(mzbwadden$ENSy, 
        main = "Effective number of species y")
hist(mzbwadden$ENSy)
dotchart(mzbwadden$ENSy)
boxplot(mzbwadden$ENS, 
        main = "Effective number of species a")
hist(mzbwadden$ENS)
dotchart(mzbwadden$ENS)
boxplot(mzbwadden$biomass, 
        main = "biomass")
hist(mzbwadden$biomass)
dotchart(mzbwadden$biomass)
par(mfrow = c(1, 1))
plot(x = mzbwadden$year, 
     y = mzbwadden$Sy)
plot(x = mzbwadden$year, 
     y = mzbwadden$ENSy)
plot(x = mzbwadden$year, 
     y = mzbwadden$S)
plot(x = mzbwadden$year, 
     y = mzbwadden$ENS)
plot(x = mzbwadden$year, 
     y = mzbwadden$biomass)

# Collinearity X
mzbwadden$biomass <- log(mzbwadden$biomass)
mzbwadden$ENS <- log(mzbwadden$ENS)
mzbwadden$ENSy <- log(mzbwadden$ENSy)
names(mzbwadden)
pairs.panels(mzbwadden[,3:10])

# Zero inflation
sum(mzbwadden == 0)

# Scale variables (method for GLMs)
mzbwadden$salinity <- scale(mzbwadden$salinity)
mzbwadden$temperature <- scale(mzbwadden$temperature)
mzbwadden$grain.size <- scale(mzbwadden$grain.size)

# Recode vars to roughly same scale
#mzbwadden$biomass <- mzbwadden$biomass/10000
mzbwadden$Sy <- mzbwadden$Sy/10
mzbwadden$S <- mzbwadden$S/10

################################# gls fitting 


gls.Sy <- gls(Sy ~ temperature + salinity + grain.size, data=mzbwadden, 
             correlation=corGaus(form= ~ lat + long | year))
summary(gls.Sy)
plot(gls.Sy)
E<-resid(gls.Sy,type="normalized")
Fit=fitted(gls.Sy)

op<-par(mfrow=c(1,2))
plot(x=Fit,y=E,xlab="Fitted values",ylab="Residuals",
     main="Residuals versus fitted values")
#identify(Fit,E)
hist(E,nclass=15)


gls.ENSy <- gls(ENSy ~ temperature + salinity, data=mzbwadden, 
               correlation=corGaus(form= ~ long + lat | year))
summary(gls.ENSy)
plot(gls.ENSy)
E<-resid(gls.ENSy,type="normalized")
Fit=fitted(gls.ENSy)

op<-par(mfrow=c(1,2))
plot(x=Fit,y=E,xlab="Fitted values",ylab="Residuals",
     main="Residuals versus fitted values")
#identify(Fit,E)
hist(E,nclass=15)


gls.biomy <- gls(biomass ~ temperature + salinity + grain.size + Sy, data=mzbwadden, 
                correlation=corGaus(form= ~ long + lat | year))
summary(gls.biomy)
plot(gls.biomy)
E<-resid(gls.biomy,type="normalized")
Fit=fitted(gls.biomy)

op<-par(mfrow=c(1,2))
plot(x=Fit,y=E,xlab="Fitted values",ylab="Residuals",
     main="Residuals versus fitted values")
#identify(Fit,E)
hist(E,nclass=15)


gls.S <- gls(S ~ temperature + salinity + grain.size, data=mzbwadden, 
              correlation=corGaus(form= ~ lat + long | year))
summary(gls.S)
plot(gls.S)
E<-resid(gls.S,type="normalized")
Fit=fitted(gls.S)
op<-par(mfrow=c(1,2))
plot(x=Fit,y=E,xlab="Fitted values",ylab="Residuals",
     main="Residuals versus fitted values")
#identify(Fit,E)
hist(E,nclass=15)


gls.ENS <- gls(ENS ~ temperature + salinity, data=mzbwadden, 
                correlation=corGaus(form= ~ long + lat | year))
summary(gls.ENS)
plot(gls.ENS)
E<-resid(gls.ENS,type="normalized")
Fit=fitted(gls.ENS)

op<-par(mfrow=c(1,2))
plot(x=Fit,y=E,xlab="Fitted values",ylab="Residuals",
     main="Residuals versus fitted values")
#identify(Fit,E)
hist(E,nclass=15)

gls.biom <- gls(biomass ~ temperature + salinity + grain.size + S, data=mzbwadden, 
                 correlation=corGaus(form= ~ long + lat | year))
summary(gls.biom)
plot(gls.biom)
E<-resid(gls.biom,type="normalized")
Fit=fitted(gls.biom)

op<-par(mfrow=c(1,2))
plot(x=Fit,y=E,xlab="Fitted values",ylab="Residuals",
     main="Residuals versus fitted values")
#identify(Fit,E)
hist(E,nclass=15)

################################# piecewiseSEM fitting

sem.SD.gls1 <- psem(gls.Sy,
                gls.ENSy,
                gls.biomy,
                Sy %~~% ENSy,
                mzbwadden)
sem.SD.gls1
summary(sem.SD.gls1)

sem.SD.gls2 <- psem(gls.S,
                   gls.ENS,
                   gls.biom,
                   S %~~% ENS,
                   mzbwadden)
sem.SD.gls2
summary(sem.SD.gls2)


####################################################### GDM Analyis

source("turnover.R")

###############################################

data <- read.csv("Leyb_MZB_3_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 7:N;


#richness based turnover index  
SERr = turnover(data[, SpecColumns], method = "SERr") #explicit
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)
new_row5 <- c(5, 5, 0)
new_row6 <- c(6, 6, 0)
new_row7 <- c(7, 7, 0)
new_row8 <- c(8, 8, 0) 
new_row9 <- c(9, 9, 0) 
SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9)
SERr2 <- acast(SERr, To ~ From , value.var="SER", na.rm= TRUE)


#abundance based turnover index 
SERa = turnover_s(data[, SpecColumns], method = "SERa")
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)
new_row5 <- c(5, 5, 0)
new_row6 <- c(6, 6, 0)
new_row7 <- c(7, 7, 0)
new_row8 <- c(8, 8, 0) 
new_row9 <- c(9, 9, 0) 
SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","grain.size","salinity"))
pred <- names(data[,env])
d <- data

################################# SERr

dcomm <- SERr2

# Dissimilarity matrix
dcomm <- cbind(year = paste(d$year), as.data.frame(dcomm))  # Adding site column for gdm function

dim(d)
dim(dcomm)


### Model ----------------------------------------------------------------------------------

# Data gdm format
print("formatting gdm data... ")
print(Sys.time())
gdm_data <- formatsitepair(bioData = dcomm, bioFormat=3, 
                           # bioData = d[,c(1,sp)], bioFormat=1, dist="bray", abundance=TRUE,
                           siteColumn="year", XColumn="long", YColumn="lat", 
                           predData = d[,c("year","long","lat", pred)])

gdm.1 <- gdm(data = gdm_data, geo = FALSE, splines = NULL, knots = NULL)

summary(gdm.1)
plot(gdm.1, plot.layout=c(2,3))

################################# SERa

dcomm <- SERa2

# Dissimilarity matrix
dcomm <- cbind(year = paste(d$year), as.data.frame(dcomm))  # Adding site column for gdm function

dim(d)
dim(dcomm)


### Model ----------------------------------------------------------------------------------

# Data gdm format
print("formatting gdm data... ")
print(Sys.time())
gdm_data <- formatsitepair(bioData = dcomm, bioFormat=3, 
                           # bioData = d[,c(1,sp)], bioFormat=1, dist="bray", abundance=TRUE,
                           siteColumn="year", XColumn="long", YColumn="lat", 
                           predData = d[,c("year","long","lat", pred)])

gdm.1 <- gdm(data = gdm_data, geo = FALSE, splines = NULL, knots = NULL)

summary(gdm.1)
plot(gdm.1, plot.layout=c(2,3))

#####################################

data <- read.csv("Nney_MZB_2_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 7:N;


#richness based turnover index  
SERr = turnover(data[, SpecColumns], method = "SERr") #explicit
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)
new_row5 <- c(5, 5, 0)
new_row6 <- c(6, 6, 0)
new_row7 <- c(7, 7, 0)
new_row8 <- c(8, 8, 0) 
new_row9 <- c(9, 9, 0) 
new_row10 <- c(10, 10, 0)
new_row11 <- c(11, 11, 0)
new_row12 <- c(12, 12, 0)
SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10, new_row11, new_row12)
SERr2 <- acast(SERr, To ~ From , value.var="SER", na.rm= TRUE)


#abundance based turnover index 
SERa = turnover_s(data[, SpecColumns], method = "SERa")
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)
new_row5 <- c(5, 5, 0)
new_row6 <- c(6, 6, 0)
new_row7 <- c(7, 7, 0)
new_row8 <- c(8, 8, 0) 
new_row9 <- c(9, 9, 0) 
new_row10 <- c(10, 10, 0)
new_row11 <- c(11, 11, 0)
new_row12 <- c(12, 12, 0)
SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10, new_row11, new_row12)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","grain.size","salinity"))
pred <- names(data[,env])
d <- data

################################# SERr

dcomm <- SERr2

# Dissimilarity matrix
dcomm <- cbind(year = paste(d$year), as.data.frame(dcomm))  # Adding site column for gdm function

dim(d)
dim(dcomm)


### Model ----------------------------------------------------------------------------------

# Data gdm format
print("formatting gdm data... ")
print(Sys.time())
gdm_data <- formatsitepair(bioData = dcomm, bioFormat=3, 
                           # bioData = d[,c(1,sp)], bioFormat=1, dist="bray", abundance=TRUE,
                           siteColumn="year", XColumn="long", YColumn="lat", 
                           predData = d[,c("year","long","lat", pred)])

gdm.1 <- gdm(data = gdm_data, geo = FALSE, splines = NULL, knots = NULL)

summary(gdm.1)
plot(gdm.1, plot.layout=c(2,3))

################################# SERa

dcomm <- SERa2

# Dissimilarity matrix
dcomm <- cbind(year = paste(d$year), as.data.frame(dcomm))  # Adding site column for gdm function

dim(d)
dim(dcomm)


### Model ----------------------------------------------------------------------------------

# Data gdm format
print("formatting gdm data... ")
print(Sys.time())
gdm_data <- formatsitepair(bioData = dcomm, bioFormat=3, 
                           # bioData = d[,c(1,sp)], bioFormat=1, dist="bray", abundance=TRUE,
                           siteColumn="year", XColumn="long", YColumn="lat", 
                           predData = d[,c("year","long","lat", pred)])

gdm.1 <- gdm(data = gdm_data, geo = FALSE, splines = NULL, knots = NULL)

summary(gdm.1)
plot(gdm.1, plot.layout=c(2,3))

#####################################

data <- read.csv("Nney_MZB_1_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 7:N;


#richness based turnover index  
SERr = turnover(data[, SpecColumns], method = "SERr") #explicit
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)
new_row5 <- c(5, 5, 0)
new_row6 <- c(6, 6, 0)
new_row7 <- c(7, 7, 0)
new_row8 <- c(8, 8, 0) 
new_row9 <- c(9, 9, 0) 
new_row10 <- c(10, 10, 0)
new_row11 <- c(11, 11, 0)
new_row12 <- c(12, 12, 0)
SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10, new_row11, new_row12)
SERr2 <- acast(SERr, To ~ From , value.var="SER", na.rm= TRUE)


#abundance based turnover index 
SERa = turnover_s(data[, SpecColumns], method = "SERa")
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)
new_row5 <- c(5, 5, 0)
new_row6 <- c(6, 6, 0)
new_row7 <- c(7, 7, 0)
new_row8 <- c(8, 8, 0) 
new_row9 <- c(9, 9, 0) 
new_row10 <- c(10, 10, 0)
new_row11 <- c(11, 11, 0)
new_row12 <- c(12, 12, 0)
SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10, new_row11, new_row12)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","grain.size","salinity"))
pred <- names(data[,env])
d <- data

################################# SERr

dcomm <- SERr2

# Dissimilarity matrix
dcomm <- cbind(year = paste(d$year), as.data.frame(dcomm))  # Adding site column for gdm function

dim(d)
dim(dcomm)


### Model ----------------------------------------------------------------------------------

# Data gdm format
print("formatting gdm data... ")
print(Sys.time())
gdm_data <- formatsitepair(bioData = dcomm, bioFormat=3, 
                           # bioData = d[,c(1,sp)], bioFormat=1, dist="bray", abundance=TRUE,
                           siteColumn="year", XColumn="long", YColumn="lat", 
                           predData = d[,c("year","long","lat", pred)])

gdm.1 <- gdm(data = gdm_data, geo = FALSE, splines = NULL, knots = NULL)

summary(gdm.1)
plot(gdm.1, plot.layout=c(2,3))

################################# SERa

dcomm <- SERa2

# Dissimilarity matrix
dcomm <- cbind(year = paste(d$year), as.data.frame(dcomm))  # Adding site column for gdm function

dim(d)
dim(dcomm)


### Model ----------------------------------------------------------------------------------

# Data gdm format
print("formatting gdm data... ")
print(Sys.time())
gdm_data <- formatsitepair(bioData = dcomm, bioFormat=3, 
                           # bioData = d[,c(1,sp)], bioFormat=1, dist="bray", abundance=TRUE,
                           siteColumn="year", XColumn="long", YColumn="lat", 
                           predData = d[,c("year","long","lat", pred)])

gdm.1 <- gdm(data = gdm_data, geo = FALSE, splines = NULL, knots = NULL)

summary(gdm.1)
plot(gdm.1, plot.layout=c(2,3))

#####################################

data <- read.csv("WeMu_MZB_3_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 7:N;

#richness based turnover index  
SERr = turnover(data[, SpecColumns], method = "SERr") #explicit
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)
new_row5 <- c(5, 5, 0)
new_row6 <- c(6, 6, 0)
new_row7 <- c(7, 7, 0)
SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7)
SERr2 <- acast(SERr, To ~ From , value.var="SER", na.rm= TRUE)


#abundance based turnover index 
SERa = turnover_s(data[, SpecColumns], method = "SERa")
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)
new_row5 <- c(5, 5, 0)
new_row6 <- c(6, 6, 0)
new_row7 <- c(7, 7, 0)
SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","grain.size","salinity"))
pred <- names(data[,env])
d <- data

################################# SERr

dcomm <- SERr2

# Dissimilarity matrix
dcomm <- cbind(year = paste(d$year), as.data.frame(dcomm))  # Adding site column for gdm function

dim(d)
dim(dcomm)


### Model ----------------------------------------------------------------------------------

# Data gdm format
print("formatting gdm data... ")
print(Sys.time())
gdm_data <- formatsitepair(bioData = dcomm, bioFormat=3, 
                           # bioData = d[,c(1,sp)], bioFormat=1, dist="bray", abundance=TRUE,
                           siteColumn="year", XColumn="long", YColumn="lat", 
                           predData = d[,c("year","long","lat", pred)])

gdm.1 <- gdm(data = gdm_data, geo = FALSE, splines = NULL, knots = NULL)

summary(gdm.1)
plot(gdm.1, plot.layout=c(2,3))

################################# SERa

dcomm <- SERa2

# Dissimilarity matrix
dcomm <- cbind(year = paste(d$year), as.data.frame(dcomm))  # Adding site column for gdm function

dim(d)
dim(dcomm)


### Model ----------------------------------------------------------------------------------

# Data gdm format
print("formatting gdm data... ")
print(Sys.time())
gdm_data <- formatsitepair(bioData = dcomm, bioFormat=3, 
                           # bioData = d[,c(1,sp)], bioFormat=1, dist="bray", abundance=TRUE,
                           siteColumn="year", XColumn="long", YColumn="lat", 
                           predData = d[,c("year","long","lat", pred)])

gdm.1 <- gdm(data = gdm_data, geo = FALSE, splines = NULL, knots = NULL)

summary(gdm.1)
plot(gdm.1, plot.layout=c(2,3))

#####################################

data <- read.csv("Nney_MZB_8_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 7:N;

#richness based turnover index  
SERr = turnover(data[, SpecColumns], method = "SERr") #explicit
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
SERr <- rbind(SERr, new_row1, new_row2, new_row3)
SERr2 <- acast(SERr, To ~ From , value.var="SER", na.rm= TRUE)


#abundance based turnover index 
SERa = turnover_s(data[, SpecColumns], method = "SERa")
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
SERa <- rbind(SERa, new_row1, new_row2, new_row3)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","grain.size","salinity"))
pred <- names(data[,env])
d <- data

################################# SERr

dcomm <- SERr2

# Dissimilarity matrix
dcomm <- cbind(year = paste(d$year), as.data.frame(dcomm))  # Adding site column for gdm function

dim(d)
dim(dcomm)


### Model ----------------------------------------------------------------------------------

# Data gdm format
print("formatting gdm data... ")
print(Sys.time())
gdm_data <- formatsitepair(bioData = dcomm, bioFormat=3, 
                           # bioData = d[,c(1,sp)], bioFormat=1, dist="bray", abundance=TRUE,
                           siteColumn="year", XColumn="long", YColumn="lat", 
                           predData = d[,c("year","long","lat", pred)])

gdm.1 <- gdm(data = gdm_data, geo = FALSE, splines = NULL, knots = NULL)

summary(gdm.1)
plot(gdm.1, plot.layout=c(2,3))

################################# SERa

dcomm <- SERa2

# Dissimilarity matrix
dcomm <- cbind(year = paste(d$year), as.data.frame(dcomm))  # Adding site column for gdm function

dim(d)
dim(dcomm)


### Model ----------------------------------------------------------------------------------

# Data gdm format
print("formatting gdm data... ")
print(Sys.time())
gdm_data <- formatsitepair(bioData = dcomm, bioFormat=3, 
                           # bioData = d[,c(1,sp)], bioFormat=1, dist="bray", abundance=TRUE,
                           siteColumn="year", XColumn="long", YColumn="lat", 
                           predData = d[,c("year","long","lat", pred)])

gdm.1 <- gdm(data = gdm_data, geo = FALSE, splines = NULL, knots = NULL)

summary(gdm.1)
plot(gdm.1, plot.layout=c(2,3))

#####################################

data <- read.csv("Nney_MZB_7_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 7:N;

#richness based turnover index  
SERr = turnover(data[, SpecColumns], method = "SERr") #explicit
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
SERr <- rbind(SERr, new_row1, new_row2, new_row3)
SERr2 <- acast(SERr, To ~ From , value.var="SER", na.rm= TRUE)


#abundance based turnover index 
SERa = turnover_s(data[, SpecColumns], method = "SERa")
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
SERa <- rbind(SERa, new_row1, new_row2, new_row3)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","grain.size","salinity"))
pred <- names(data[,env])
d <- data

################################# SERr

dcomm <- SERr2

# Dissimilarity matrix
dcomm <- cbind(year = paste(d$year), as.data.frame(dcomm))  # Adding site column for gdm function

dim(d)
dim(dcomm)


### Model ----------------------------------------------------------------------------------

# Data gdm format
print("formatting gdm data... ")
print(Sys.time())
gdm_data <- formatsitepair(bioData = dcomm, bioFormat=3, 
                           # bioData = d[,c(1,sp)], bioFormat=1, dist="bray", abundance=TRUE,
                           siteColumn="year", XColumn="long", YColumn="lat", 
                           predData = d[,c("year","long","lat", pred)])

gdm.1 <- gdm(data = gdm_data, geo = FALSE, splines = NULL, knots = NULL)

summary(gdm.1)
plot(gdm.1, plot.layout=c(2,3))

################################# SERa

dcomm <- SERa2

# Dissimilarity matrix
dcomm <- cbind(year = paste(d$year), as.data.frame(dcomm))  # Adding site column for gdm function

dim(d)
dim(dcomm)


### Model ----------------------------------------------------------------------------------

# Data gdm format
print("formatting gdm data... ")
print(Sys.time())
gdm_data <- formatsitepair(bioData = dcomm, bioFormat=3, 
                           # bioData = d[,c(1,sp)], bioFormat=1, dist="bray", abundance=TRUE,
                           siteColumn="year", XColumn="long", YColumn="lat", 
                           predData = d[,c("year","long","lat", pred)])

gdm.1 <- gdm(data = gdm_data, geo = FALSE, splines = NULL, knots = NULL)

summary(gdm.1)
plot(gdm.1, plot.layout=c(2,3))

#####################################

data <- read.csv("Nney_MZB_6_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 7:N;

#richness based turnover index  
SERr = turnover(data[, SpecColumns], method = "SERr") #explicit
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
SERr <- rbind(SERr, new_row1, new_row2, new_row3)
SERr2 <- acast(SERr, To ~ From , value.var="SER", na.rm= TRUE)


#abundance based turnover index 
SERa = turnover_s(data[, SpecColumns], method = "SERa")
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
SERa <- rbind(SERa, new_row1, new_row2, new_row3)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","grain.size","salinity"))
pred <- names(data[,env])
d <- data

################################# SERr

dcomm <- SERr2

# Dissimilarity matrix
dcomm <- cbind(year = paste(d$year), as.data.frame(dcomm))  # Adding site column for gdm function

dim(d)
dim(dcomm)


### Model ----------------------------------------------------------------------------------

# Data gdm format
print("formatting gdm data... ")
print(Sys.time())
gdm_data <- formatsitepair(bioData = dcomm, bioFormat=3, 
                           # bioData = d[,c(1,sp)], bioFormat=1, dist="bray", abundance=TRUE,
                           siteColumn="year", XColumn="long", YColumn="lat", 
                           predData = d[,c("year","long","lat", pred)])

gdm.1 <- gdm(data = gdm_data, geo = FALSE, splines = NULL, knots = NULL)

summary(gdm.1)
plot(gdm.1, plot.layout=c(2,3))

################################# SERa

dcomm <- SERa2

# Dissimilarity matrix
dcomm <- cbind(year = paste(d$year), as.data.frame(dcomm))  # Adding site column for gdm function

dim(d)
dim(dcomm)


### Model ----------------------------------------------------------------------------------

# Data gdm format
print("formatting gdm data... ")
print(Sys.time())
gdm_data <- formatsitepair(bioData = dcomm, bioFormat=3, 
                           # bioData = d[,c(1,sp)], bioFormat=1, dist="bray", abundance=TRUE,
                           siteColumn="year", XColumn="long", YColumn="lat", 
                           predData = d[,c("year","long","lat", pred)])

gdm.1 <- gdm(data = gdm_data, geo = FALSE, splines = NULL, knots = NULL)

summary(gdm.1)
plot(gdm.1, plot.layout=c(2,3))

#####################################

data <- read.csv("Nney_MZB_3_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 7:N;


#richness based turnover index  
SERr = turnover(data[, SpecColumns], method = "SERr") #explicit
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)
new_row5 <- c(5, 5, 0)
new_row6 <- c(6, 6, 0)
new_row7 <- c(7, 7, 0)
new_row8 <- c(8, 8, 0) 
new_row9 <- c(9, 9, 0) 
new_row10 <- c(10, 10, 0)
new_row11 <- c(11, 11, 0)
SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10, new_row11)
SERr2 <- acast(SERr, To ~ From , value.var="SER", na.rm= TRUE)


#abundance based turnover index 
SERa = turnover_s(data[, SpecColumns], method = "SERa")
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)
new_row5 <- c(5, 5, 0)
new_row6 <- c(6, 6, 0)
new_row7 <- c(7, 7, 0)
new_row8 <- c(8, 8, 0) 
new_row9 <- c(9, 9, 0) 
new_row10 <- c(10, 10, 0)
new_row11 <- c(11, 11, 0)
SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10, new_row11)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","grain.size","salinity"))
pred <- names(data[,env])
d <- data

################################# SERr

dcomm <- SERr2

# Dissimilarity matrix
dcomm <- cbind(year = paste(d$year), as.data.frame(dcomm))  # Adding site column for gdm function

dim(d)
dim(dcomm)


### Model ----------------------------------------------------------------------------------

# Data gdm format
print("formatting gdm data... ")
print(Sys.time())
gdm_data <- formatsitepair(bioData = dcomm, bioFormat=3, 
                           # bioData = d[,c(1,sp)], bioFormat=1, dist="bray", abundance=TRUE,
                           siteColumn="year", XColumn="long", YColumn="lat", 
                           predData = d[,c("year","long","lat", pred)])

gdm.1 <- gdm(data = gdm_data, geo = FALSE, splines = NULL, knots = NULL)

summary(gdm.1)
plot(gdm.1, plot.layout=c(2,3))

################################# SERa

dcomm <- SERa2

# Dissimilarity matrix
dcomm <- cbind(year = paste(d$year), as.data.frame(dcomm))  # Adding site column for gdm function

dim(d)
dim(dcomm)


### Model ----------------------------------------------------------------------------------

# Data gdm format
print("formatting gdm data... ")
print(Sys.time())
gdm_data <- formatsitepair(bioData = dcomm, bioFormat=3, 
                           # bioData = d[,c(1,sp)], bioFormat=1, dist="bray", abundance=TRUE,
                           siteColumn="year", XColumn="long", YColumn="lat", 
                           predData = d[,c("year","long","lat", pred)])

gdm.1 <- gdm(data = gdm_data, geo = FALSE, splines = NULL, knots = NULL)

summary(gdm.1)
plot(gdm.1, plot.layout=c(2,3))


#####################################

data <- read.csv("1_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 7:N;


#richness based turnover index  
SERr = turnover(data[, SpecColumns], method = "SERr") #explicit
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)
new_row5 <- c(5, 5, 0)
new_row6 <- c(6, 6, 0)
new_row7 <- c(7, 7, 0)
new_row8 <- c(8, 8, 0) 
new_row9 <- c(9, 9, 0) 
new_row10 <- c(10, 10, 0)
new_row11 <- c(11, 11, 0)
SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10, new_row11)
SERr2 <- acast(SERr, To ~ From , value.var="SER", na.rm= TRUE)


#abundance based turnover index 
SERa = turnover_s(data[, SpecColumns], method = "SERa")
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)
new_row5 <- c(5, 5, 0)
new_row6 <- c(6, 6, 0)
new_row7 <- c(7, 7, 0)
new_row8 <- c(8, 8, 0) 
new_row9 <- c(9, 9, 0) 
new_row10 <- c(10, 10, 0)
new_row11 <- c(11, 11, 0)
SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10, new_row11)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","grain.size","salinity"))
pred <- names(data[,env])
d <- data

################################# SERr

dcomm <- SERr2

# Dissimilarity matrix
dcomm <- cbind(year = paste(d$year), as.data.frame(dcomm))  # Adding site column for gdm function

dim(d)
dim(dcomm)


### Model ----------------------------------------------------------------------------------

# Data gdm format
print("formatting gdm data... ")
print(Sys.time())
gdm_data <- formatsitepair(bioData = dcomm, bioFormat=3, 
                           # bioData = d[,c(1,sp)], bioFormat=1, dist="bray", abundance=TRUE,
                           siteColumn="year", XColumn="long", YColumn="lat", 
                           predData = d[,c("year","long","lat", pred)])

gdm.1 <- gdm(data = gdm_data, geo = FALSE, splines = NULL, knots = NULL)

summary(gdm.1)
plot(gdm.1, plot.layout=c(2,3))

################################# SERa

dcomm <- SERa2

# Dissimilarity matrix
dcomm <- cbind(year = paste(d$year), as.data.frame(dcomm))  # Adding site column for gdm function

dim(d)
dim(dcomm)


### Model ----------------------------------------------------------------------------------

# Data gdm format
print("formatting gdm data... ")
print(Sys.time())
gdm_data <- formatsitepair(bioData = dcomm, bioFormat=3, 
                           # bioData = d[,c(1,sp)], bioFormat=1, dist="bray", abundance=TRUE,
                           siteColumn="year", XColumn="long", YColumn="lat", 
                           predData = d[,c("year","long","lat", pred)])

gdm.1 <- gdm(data = gdm_data, geo = FALSE, splines = NULL, knots = NULL)

summary(gdm.1)
plot(gdm.1, plot.layout=c(2,3))

#####################################

data <- read.csv("10_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 7:N;


#richness based turnover index  
SERr = turnover(data[, SpecColumns], method = "SERr") #explicit
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)
new_row5 <- c(5, 5, 0)
new_row6 <- c(6, 6, 0)
new_row7 <- c(7, 7, 0)
new_row8 <- c(8, 8, 0) 
new_row9 <- c(9, 9, 0) 
new_row10 <- c(10, 10, 0)
new_row11 <- c(11, 11, 0)
SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10, new_row11)
SERr2 <- acast(SERr, To ~ From , value.var="SER", na.rm= TRUE)


#abundance based turnover index 
SERa = turnover_s(data[, SpecColumns], method = "SERa")
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)
new_row5 <- c(5, 5, 0)
new_row6 <- c(6, 6, 0)
new_row7 <- c(7, 7, 0)
new_row8 <- c(8, 8, 0) 
new_row9 <- c(9, 9, 0) 
new_row10 <- c(10, 10, 0)
new_row11 <- c(11, 11, 0)
SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10, new_row11)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","grain.size","salinity"))
pred <- names(data[,env])
d <- data

################################# SERr

dcomm <- SERr2

# Dissimilarity matrix
dcomm <- cbind(year = paste(d$year), as.data.frame(dcomm))  # Adding site column for gdm function

dim(d)
dim(dcomm)


### Model ----------------------------------------------------------------------------------

# Data gdm format
print("formatting gdm data... ")
print(Sys.time())
gdm_data <- formatsitepair(bioData = dcomm, bioFormat=3, 
                           # bioData = d[,c(1,sp)], bioFormat=1, dist="bray", abundance=TRUE,
                           siteColumn="year", XColumn="long", YColumn="lat", 
                           predData = d[,c("year","long","lat", pred)])

gdm.1 <- gdm(data = gdm_data, geo = FALSE, splines = NULL, knots = NULL)

summary(gdm.1)
plot(gdm.1, plot.layout=c(2,3))

################################# SERa

dcomm <- SERa2

# Dissimilarity matrix
dcomm <- cbind(year = paste(d$year), as.data.frame(dcomm))  # Adding site column for gdm function

dim(d)
dim(dcomm)


### Model ----------------------------------------------------------------------------------

# Data gdm format
print("formatting gdm data... ")
print(Sys.time())
gdm_data <- formatsitepair(bioData = dcomm, bioFormat=3, 
                           # bioData = d[,c(1,sp)], bioFormat=1, dist="bray", abundance=TRUE,
                           siteColumn="year", XColumn="long", YColumn="lat", 
                           predData = d[,c("year","long","lat", pred)])

gdm.1 <- gdm(data = gdm_data, geo = FALSE, splines = NULL, knots = NULL)

summary(gdm.1)
plot(gdm.1, plot.layout=c(2,3))

#####################################

data <- read.csv("11_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 7:N;


#richness based turnover index  
SERr = turnover(data[, SpecColumns], method = "SERr") #explicit
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)
new_row5 <- c(5, 5, 0)
new_row6 <- c(6, 6, 0)
new_row7 <- c(7, 7, 0)
new_row8 <- c(8, 8, 0) 
new_row9 <- c(9, 9, 0) 
new_row10 <- c(10, 10, 0)
new_row11 <- c(11, 11, 0)
SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10, new_row11)
SERr2 <- acast(SERr, To ~ From , value.var="SER", na.rm= TRUE)


#abundance based turnover index 
SERa = turnover_s(data[, SpecColumns], method = "SERa")
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)
new_row5 <- c(5, 5, 0)
new_row6 <- c(6, 6, 0)
new_row7 <- c(7, 7, 0)
new_row8 <- c(8, 8, 0) 
new_row9 <- c(9, 9, 0) 
new_row10 <- c(10, 10, 0)
new_row11 <- c(11, 11, 0)
SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10, new_row11)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","grain.size","salinity"))
pred <- names(data[,env])
d <- data

################################# SERr

dcomm <- SERr2

# Dissimilarity matrix
dcomm <- cbind(year = paste(d$year), as.data.frame(dcomm))  # Adding site column for gdm function

dim(d)
dim(dcomm)


### Model ----------------------------------------------------------------------------------

# Data gdm format
print("formatting gdm data... ")
print(Sys.time())
gdm_data <- formatsitepair(bioData = dcomm, bioFormat=3, 
                           # bioData = d[,c(1,sp)], bioFormat=1, dist="bray", abundance=TRUE,
                           siteColumn="year", XColumn="long", YColumn="lat", 
                           predData = d[,c("year","long","lat", pred)])

gdm.1 <- gdm(data = gdm_data, geo = FALSE, splines = NULL, knots = NULL)

summary(gdm.1)
plot(gdm.1, plot.layout=c(2,3))

################################# SERa

dcomm <- SERa2

# Dissimilarity matrix
dcomm <- cbind(year = paste(d$year), as.data.frame(dcomm))  # Adding site column for gdm function

dim(d)
dim(dcomm)


### Model ----------------------------------------------------------------------------------

# Data gdm format
print("formatting gdm data... ")
print(Sys.time())
gdm_data <- formatsitepair(bioData = dcomm, bioFormat=3, 
                           # bioData = d[,c(1,sp)], bioFormat=1, dist="bray", abundance=TRUE,
                           siteColumn="year", XColumn="long", YColumn="lat", 
                           predData = d[,c("year","long","lat", pred)])

gdm.1 <- gdm(data = gdm_data, geo = FALSE, splines = NULL, knots = NULL)

summary(gdm.1)
plot(gdm.1, plot.layout=c(2,3))

#####################################

data <- read.csv("12_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 7:N;


#richness based turnover index  
SERr = turnover(data[, SpecColumns], method = "SERr") #explicit
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)
new_row5 <- c(5, 5, 0)
new_row6 <- c(6, 6, 0)
new_row7 <- c(7, 7, 0)
new_row8 <- c(8, 8, 0) 
new_row9 <- c(9, 9, 0) 
new_row10 <- c(10, 10, 0)
new_row11 <- c(11, 11, 0)
SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10, new_row11)
SERr2 <- acast(SERr, To ~ From , value.var="SER", na.rm= TRUE)


#abundance based turnover index 
SERa = turnover_s(data[, SpecColumns], method = "SERa")
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)
new_row5 <- c(5, 5, 0)
new_row6 <- c(6, 6, 0)
new_row7 <- c(7, 7, 0)
new_row8 <- c(8, 8, 0) 
new_row9 <- c(9, 9, 0) 
new_row10 <- c(10, 10, 0)
new_row11 <- c(11, 11, 0)
SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10, new_row11)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","grain.size","salinity"))
pred <- names(data[,env])
d <- data

################################# SERr

dcomm <- SERr2

# Dissimilarity matrix
dcomm <- cbind(year = paste(d$year), as.data.frame(dcomm))  # Adding site column for gdm function

dim(d)
dim(dcomm)


### Model ----------------------------------------------------------------------------------

# Data gdm format
print("formatting gdm data... ")
print(Sys.time())
gdm_data <- formatsitepair(bioData = dcomm, bioFormat=3, 
                           # bioData = d[,c(1,sp)], bioFormat=1, dist="bray", abundance=TRUE,
                           siteColumn="year", XColumn="long", YColumn="lat", 
                           predData = d[,c("year","long","lat", pred)])

gdm.1 <- gdm(data = gdm_data, geo = FALSE, splines = NULL, knots = NULL)

summary(gdm.1)
plot(gdm.1, plot.layout=c(2,3))

################################# SERa

dcomm <- SERa2

# Dissimilarity matrix
dcomm <- cbind(year = paste(d$year), as.data.frame(dcomm))  # Adding site column for gdm function

dim(d)
dim(dcomm)


### Model ----------------------------------------------------------------------------------

# Data gdm format
print("formatting gdm data... ")
print(Sys.time())
gdm_data <- formatsitepair(bioData = dcomm, bioFormat=3, 
                           # bioData = d[,c(1,sp)], bioFormat=1, dist="bray", abundance=TRUE,
                           siteColumn="year", XColumn="long", YColumn="lat", 
                           predData = d[,c("year","long","lat", pred)])

gdm.1 <- gdm(data = gdm_data, geo = FALSE, splines = NULL, knots = NULL)

summary(gdm.1)
plot(gdm.1, plot.layout=c(2,3))

#####################################

data <- read.csv("2_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 7:N;


#richness based turnover index  
SERr = turnover(data[, SpecColumns], method = "SERr") #explicit
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)
new_row5 <- c(5, 5, 0)
new_row6 <- c(6, 6, 0)
new_row7 <- c(7, 7, 0)
new_row8 <- c(8, 8, 0) 
new_row9 <- c(9, 9, 0) 
new_row10 <- c(10, 10, 0)
new_row11 <- c(11, 11, 0)
SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10, new_row11)
SERr2 <- acast(SERr, To ~ From , value.var="SER", na.rm= TRUE)


#abundance based turnover index 
SERa = turnover_s(data[, SpecColumns], method = "SERa")
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)
new_row5 <- c(5, 5, 0)
new_row6 <- c(6, 6, 0)
new_row7 <- c(7, 7, 0)
new_row8 <- c(8, 8, 0) 
new_row9 <- c(9, 9, 0) 
new_row10 <- c(10, 10, 0)
new_row11 <- c(11, 11, 0)
SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10, new_row11)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","grain.size","salinity"))
pred <- names(data[,env])
d <- data

################################# SERr

dcomm <- SERr2

# Dissimilarity matrix
dcomm <- cbind(year = paste(d$year), as.data.frame(dcomm))  # Adding site column for gdm function

dim(d)
dim(dcomm)


### Model ----------------------------------------------------------------------------------

# Data gdm format
print("formatting gdm data... ")
print(Sys.time())
gdm_data <- formatsitepair(bioData = dcomm, bioFormat=3, 
                           # bioData = d[,c(1,sp)], bioFormat=1, dist="bray", abundance=TRUE,
                           siteColumn="year", XColumn="long", YColumn="lat", 
                           predData = d[,c("year","long","lat", pred)])

gdm.1 <- gdm(data = gdm_data, geo = FALSE, splines = NULL, knots = NULL)

summary(gdm.1)
plot(gdm.1, plot.layout=c(2,3))

################################# SERa

dcomm <- SERa2

# Dissimilarity matrix
dcomm <- cbind(year = paste(d$year), as.data.frame(dcomm))  # Adding site column for gdm function

dim(d)
dim(dcomm)


### Model ----------------------------------------------------------------------------------

# Data gdm format
print("formatting gdm data... ")
print(Sys.time())
gdm_data <- formatsitepair(bioData = dcomm, bioFormat=3, 
                           # bioData = d[,c(1,sp)], bioFormat=1, dist="bray", abundance=TRUE,
                           siteColumn="year", XColumn="long", YColumn="lat", 
                           predData = d[,c("year","long","lat", pred)])

gdm.1 <- gdm(data = gdm_data, geo = FALSE, splines = NULL, knots = NULL)

summary(gdm.1)
plot(gdm.1, plot.layout=c(2,3))

#####################################

data <- read.csv("3_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 7:N;


#richness based turnover index  
SERr = turnover(data[, SpecColumns], method = "SERr") #explicit
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)
new_row5 <- c(5, 5, 0)
new_row6 <- c(6, 6, 0)
new_row7 <- c(7, 7, 0)
new_row8 <- c(8, 8, 0) 
new_row9 <- c(9, 9, 0) 
new_row10 <- c(10, 10, 0)
new_row11 <- c(11, 11, 0)
SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10, new_row11)
SERr2 <- acast(SERr, To ~ From , value.var="SER", na.rm= TRUE)


#abundance based turnover index 
SERa = turnover_s(data[, SpecColumns], method = "SERa")
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)
new_row5 <- c(5, 5, 0)
new_row6 <- c(6, 6, 0)
new_row7 <- c(7, 7, 0)
new_row8 <- c(8, 8, 0) 
new_row9 <- c(9, 9, 0) 
new_row10 <- c(10, 10, 0)
new_row11 <- c(11, 11, 0)
SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10, new_row11)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","grain.size","salinity"))
pred <- names(data[,env])
d <- data

################################# SERr

dcomm <- SERr2

# Dissimilarity matrix
dcomm <- cbind(year = paste(d$year), as.data.frame(dcomm))  # Adding site column for gdm function

dim(d)
dim(dcomm)


### Model ----------------------------------------------------------------------------------

# Data gdm format
print("formatting gdm data... ")
print(Sys.time())
gdm_data <- formatsitepair(bioData = dcomm, bioFormat=3, 
                           # bioData = d[,c(1,sp)], bioFormat=1, dist="bray", abundance=TRUE,
                           siteColumn="year", XColumn="long", YColumn="lat", 
                           predData = d[,c("year","long","lat", pred)])

gdm.1 <- gdm(data = gdm_data, geo = FALSE, splines = NULL, knots = NULL)

summary(gdm.1)
plot(gdm.1, plot.layout=c(2,3))

################################# SERa

dcomm <- SERa2

# Dissimilarity matrix
dcomm <- cbind(year = paste(d$year), as.data.frame(dcomm))  # Adding site column for gdm function

dim(d)
dim(dcomm)


### Model ----------------------------------------------------------------------------------

# Data gdm format
print("formatting gdm data... ")
print(Sys.time())
gdm_data <- formatsitepair(bioData = dcomm, bioFormat=3, 
                           # bioData = d[,c(1,sp)], bioFormat=1, dist="bray", abundance=TRUE,
                           siteColumn="year", XColumn="long", YColumn="lat", 
                           predData = d[,c("year","long","lat", pred)])

gdm.1 <- gdm(data = gdm_data, geo = FALSE, splines = NULL, knots = NULL)

summary(gdm.1)
plot(gdm.1, plot.layout=c(2,3))

#####################################

data <- read.csv("4_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 7:N;


#richness based turnover index  
SERr = turnover(data[, SpecColumns], method = "SERr") #explicit
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)
new_row5 <- c(5, 5, 0)
new_row6 <- c(6, 6, 0)
new_row7 <- c(7, 7, 0)
new_row8 <- c(8, 8, 0) 
new_row9 <- c(9, 9, 0) 
new_row10 <- c(10, 10, 0)
new_row11 <- c(11, 11, 0)
SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10, new_row11)
SERr2 <- acast(SERr, To ~ From , value.var="SER", na.rm= TRUE)


#abundance based turnover index 
SERa = turnover_s(data[, SpecColumns], method = "SERa")
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)
new_row5 <- c(5, 5, 0)
new_row6 <- c(6, 6, 0)
new_row7 <- c(7, 7, 0)
new_row8 <- c(8, 8, 0) 
new_row9 <- c(9, 9, 0) 
new_row10 <- c(10, 10, 0)
new_row11 <- c(11, 11, 0)
SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10, new_row11)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","grain.size","salinity"))
pred <- names(data[,env])
d <- data

################################# SERr

dcomm <- SERr2

# Dissimilarity matrix
dcomm <- cbind(year = paste(d$year), as.data.frame(dcomm))  # Adding site column for gdm function

dim(d)
dim(dcomm)


### Model ----------------------------------------------------------------------------------

# Data gdm format
print("formatting gdm data... ")
print(Sys.time())
gdm_data <- formatsitepair(bioData = dcomm, bioFormat=3, 
                           # bioData = d[,c(1,sp)], bioFormat=1, dist="bray", abundance=TRUE,
                           siteColumn="year", XColumn="long", YColumn="lat", 
                           predData = d[,c("year","long","lat", pred)])

gdm.1 <- gdm(data = gdm_data, geo = FALSE, splines = NULL, knots = NULL)

summary(gdm.1)
plot(gdm.1, plot.layout=c(2,3))

################################# SERa

dcomm <- SERa2

# Dissimilarity matrix
dcomm <- cbind(year = paste(d$year), as.data.frame(dcomm))  # Adding site column for gdm function

dim(d)
dim(dcomm)


### Model ----------------------------------------------------------------------------------

# Data gdm format
print("formatting gdm data... ")
print(Sys.time())
gdm_data <- formatsitepair(bioData = dcomm, bioFormat=3, 
                           # bioData = d[,c(1,sp)], bioFormat=1, dist="bray", abundance=TRUE,
                           siteColumn="year", XColumn="long", YColumn="lat", 
                           predData = d[,c("year","long","lat", pred)])

gdm.1 <- gdm(data = gdm_data, geo = FALSE, splines = NULL, knots = NULL)

summary(gdm.1)
plot(gdm.1, plot.layout=c(2,3))

#####################################

data <- read.csv("5_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 7:N;


#richness based turnover index  
SERr = turnover(data[, SpecColumns], method = "SERr") #explicit
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)
new_row5 <- c(5, 5, 0)
new_row6 <- c(6, 6, 0)
new_row7 <- c(7, 7, 0)
new_row8 <- c(8, 8, 0) 
new_row9 <- c(9, 9, 0) 
new_row10 <- c(10, 10, 0)
new_row11 <- c(11, 11, 0)
SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10, new_row11)
SERr2 <- acast(SERr, To ~ From , value.var="SER", na.rm= TRUE)


#abundance based turnover index 
SERa = turnover_s(data[, SpecColumns], method = "SERa")
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)
new_row5 <- c(5, 5, 0)
new_row6 <- c(6, 6, 0)
new_row7 <- c(7, 7, 0)
new_row8 <- c(8, 8, 0) 
new_row9 <- c(9, 9, 0) 
new_row10 <- c(10, 10, 0)
new_row11 <- c(11, 11, 0)
SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10, new_row11)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","grain.size","salinity"))
pred <- names(data[,env])
d <- data

################################# SERr

dcomm <- SERr2

# Dissimilarity matrix
dcomm <- cbind(year = paste(d$year), as.data.frame(dcomm))  # Adding site column for gdm function

dim(d)
dim(dcomm)


### Model ----------------------------------------------------------------------------------

# Data gdm format
print("formatting gdm data... ")
print(Sys.time())
gdm_data <- formatsitepair(bioData = dcomm, bioFormat=3, 
                           # bioData = d[,c(1,sp)], bioFormat=1, dist="bray", abundance=TRUE,
                           siteColumn="year", XColumn="long", YColumn="lat", 
                           predData = d[,c("year","long","lat", pred)])

gdm.1 <- gdm(data = gdm_data, geo = FALSE, splines = NULL, knots = NULL)

summary(gdm.1)
plot(gdm.1, plot.layout=c(2,3))

################################# SERa

dcomm <- SERa2

# Dissimilarity matrix
dcomm <- cbind(year = paste(d$year), as.data.frame(dcomm))  # Adding site column for gdm function

dim(d)
dim(dcomm)


### Model ----------------------------------------------------------------------------------

# Data gdm format
print("formatting gdm data... ")
print(Sys.time())
gdm_data <- formatsitepair(bioData = dcomm, bioFormat=3, 
                           # bioData = d[,c(1,sp)], bioFormat=1, dist="bray", abundance=TRUE,
                           siteColumn="year", XColumn="long", YColumn="lat", 
                           predData = d[,c("year","long","lat", pred)])

gdm.1 <- gdm(data = gdm_data, geo = FALSE, splines = NULL, knots = NULL)

summary(gdm.1)
plot(gdm.1, plot.layout=c(2,3))

#####################################

data <- read.csv("6_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 7:N;


#richness based turnover index  
SERr = turnover(data[, SpecColumns], method = "SERr") #explicit
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)
new_row5 <- c(5, 5, 0)
new_row6 <- c(6, 6, 0)
new_row7 <- c(7, 7, 0)
new_row8 <- c(8, 8, 0) 
new_row9 <- c(9, 9, 0) 
new_row10 <- c(10, 10, 0)
new_row11 <- c(11, 11, 0)
SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10, new_row11)
SERr2 <- acast(SERr, To ~ From , value.var="SER", na.rm= TRUE)


#abundance based turnover index 
SERa = turnover_s(data[, SpecColumns], method = "SERa")
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)
new_row5 <- c(5, 5, 0)
new_row6 <- c(6, 6, 0)
new_row7 <- c(7, 7, 0)
new_row8 <- c(8, 8, 0) 
new_row9 <- c(9, 9, 0) 
new_row10 <- c(10, 10, 0)
new_row11 <- c(11, 11, 0)
SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10, new_row11)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","grain.size","salinity"))
pred <- names(data[,env])
d <- data

################################# SERr

dcomm <- SERr2

# Dissimilarity matrix
dcomm <- cbind(year = paste(d$year), as.data.frame(dcomm))  # Adding site column for gdm function

dim(d)
dim(dcomm)


### Model ----------------------------------------------------------------------------------

# Data gdm format
print("formatting gdm data... ")
print(Sys.time())
gdm_data <- formatsitepair(bioData = dcomm, bioFormat=3, 
                           # bioData = d[,c(1,sp)], bioFormat=1, dist="bray", abundance=TRUE,
                           siteColumn="year", XColumn="long", YColumn="lat", 
                           predData = d[,c("year","long","lat", pred)])

gdm.1 <- gdm(data = gdm_data, geo = FALSE, splines = NULL, knots = NULL)

summary(gdm.1)
plot(gdm.1, plot.layout=c(2,3))

################################# SERa

dcomm <- SERa2

# Dissimilarity matrix
dcomm <- cbind(year = paste(d$year), as.data.frame(dcomm))  # Adding site column for gdm function

dim(d)
dim(dcomm)


### Model ----------------------------------------------------------------------------------

# Data gdm format
print("formatting gdm data... ")
print(Sys.time())
gdm_data <- formatsitepair(bioData = dcomm, bioFormat=3, 
                           # bioData = d[,c(1,sp)], bioFormat=1, dist="bray", abundance=TRUE,
                           siteColumn="year", XColumn="long", YColumn="lat", 
                           predData = d[,c("year","long","lat", pred)])

gdm.1 <- gdm(data = gdm_data, geo = FALSE, splines = NULL, knots = NULL)

summary(gdm.1)
plot(gdm.1, plot.layout=c(2,3))

#####################################

data <- read.csv("7_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 7:N;


#richness based turnover index  
SERr = turnover(data[, SpecColumns], method = "SERr") #explicit
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)
new_row5 <- c(5, 5, 0)
new_row6 <- c(6, 6, 0)
new_row7 <- c(7, 7, 0)
new_row8 <- c(8, 8, 0) 
new_row9 <- c(9, 9, 0) 
new_row10 <- c(10, 10, 0)
new_row11 <- c(11, 11, 0)
SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10, new_row11)
SERr2 <- acast(SERr, To ~ From , value.var="SER", na.rm= TRUE)


#abundance based turnover index 
SERa = turnover_s(data[, SpecColumns], method = "SERa")
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)
new_row5 <- c(5, 5, 0)
new_row6 <- c(6, 6, 0)
new_row7 <- c(7, 7, 0)
new_row8 <- c(8, 8, 0) 
new_row9 <- c(9, 9, 0) 
new_row10 <- c(10, 10, 0)
new_row11 <- c(11, 11, 0)
SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10, new_row11)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","grain.size","salinity"))
pred <- names(data[,env])
d <- data

################################# SERr

dcomm <- SERr2

# Dissimilarity matrix
dcomm <- cbind(year = paste(d$year), as.data.frame(dcomm))  # Adding site column for gdm function

dim(d)
dim(dcomm)


### Model ----------------------------------------------------------------------------------

# Data gdm format
print("formatting gdm data... ")
print(Sys.time())
gdm_data <- formatsitepair(bioData = dcomm, bioFormat=3, 
                           # bioData = d[,c(1,sp)], bioFormat=1, dist="bray", abundance=TRUE,
                           siteColumn="year", XColumn="long", YColumn="lat", 
                           predData = d[,c("year","long","lat", pred)])

gdm.1 <- gdm(data = gdm_data, geo = FALSE, splines = NULL, knots = NULL)

summary(gdm.1)
plot(gdm.1, plot.layout=c(2,3))

################################# SERa

dcomm <- SERa2

# Dissimilarity matrix
dcomm <- cbind(year = paste(d$year), as.data.frame(dcomm))  # Adding site column for gdm function

dim(d)
dim(dcomm)


### Model ----------------------------------------------------------------------------------

# Data gdm format
print("formatting gdm data... ")
print(Sys.time())
gdm_data <- formatsitepair(bioData = dcomm, bioFormat=3, 
                           # bioData = d[,c(1,sp)], bioFormat=1, dist="bray", abundance=TRUE,
                           siteColumn="year", XColumn="long", YColumn="lat", 
                           predData = d[,c("year","long","lat", pred)])

gdm.1 <- gdm(data = gdm_data, geo = FALSE, splines = NULL, knots = NULL)

summary(gdm.1)
plot(gdm.1, plot.layout=c(2,3))

#####################################

data <- read.csv("8_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 7:N;


#richness based turnover index  
SERr = turnover(data[, SpecColumns], method = "SERr") #explicit
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)
new_row5 <- c(5, 5, 0)
new_row6 <- c(6, 6, 0)
new_row7 <- c(7, 7, 0)
new_row8 <- c(8, 8, 0) 
new_row9 <- c(9, 9, 0) 
new_row10 <- c(10, 10, 0)
new_row11 <- c(11, 11, 0)
SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10, new_row11)
SERr2 <- acast(SERr, To ~ From , value.var="SER", na.rm= TRUE)


#abundance based turnover index 
SERa = turnover_s(data[, SpecColumns], method = "SERa")
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)
new_row5 <- c(5, 5, 0)
new_row6 <- c(6, 6, 0)
new_row7 <- c(7, 7, 0)
new_row8 <- c(8, 8, 0) 
new_row9 <- c(9, 9, 0) 
new_row10 <- c(10, 10, 0)
new_row11 <- c(11, 11, 0)
SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10, new_row11)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","grain.size","salinity"))
pred <- names(data[,env])
d <- data

################################# SERr

dcomm <- SERr2

# Dissimilarity matrix
dcomm <- cbind(year = paste(d$year), as.data.frame(dcomm))  # Adding site column for gdm function

dim(d)
dim(dcomm)


### Model ----------------------------------------------------------------------------------

# Data gdm format
print("formatting gdm data... ")
print(Sys.time())
gdm_data <- formatsitepair(bioData = dcomm, bioFormat=3, 
                           # bioData = d[,c(1,sp)], bioFormat=1, dist="bray", abundance=TRUE,
                           siteColumn="year", XColumn="long", YColumn="lat", 
                           predData = d[,c("year","long","lat", pred)])

gdm.1 <- gdm(data = gdm_data, geo = FALSE, splines = NULL, knots = NULL)

summary(gdm.1)
plot(gdm.1, plot.layout=c(2,3))

################################# SERa

dcomm <- SERa2

# Dissimilarity matrix
dcomm <- cbind(year = paste(d$year), as.data.frame(dcomm))  # Adding site column for gdm function

dim(d)
dim(dcomm)


### Model ----------------------------------------------------------------------------------

# Data gdm format
print("formatting gdm data... ")
print(Sys.time())
gdm_data <- formatsitepair(bioData = dcomm, bioFormat=3, 
                           # bioData = d[,c(1,sp)], bioFormat=1, dist="bray", abundance=TRUE,
                           siteColumn="year", XColumn="long", YColumn="lat", 
                           predData = d[,c("year","long","lat", pred)])

gdm.1 <- gdm(data = gdm_data, geo = FALSE, splines = NULL, knots = NULL)

summary(gdm.1)
plot(gdm.1, plot.layout=c(2,3))

#####################################

data <- read.csv("9_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 7:N;


#richness based turnover index  
SERr = turnover(data[, SpecColumns], method = "SERr") #explicit
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)
new_row5 <- c(5, 5, 0)
new_row6 <- c(6, 6, 0)
new_row7 <- c(7, 7, 0)
new_row8 <- c(8, 8, 0) 
new_row9 <- c(9, 9, 0) 
new_row10 <- c(10, 10, 0)
new_row11 <- c(11, 11, 0)
SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10, new_row11)
SERr2 <- acast(SERr, To ~ From , value.var="SER", na.rm= TRUE)


#abundance based turnover index 
SERa = turnover_s(data[, SpecColumns], method = "SERa")
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)
new_row5 <- c(5, 5, 0)
new_row6 <- c(6, 6, 0)
new_row7 <- c(7, 7, 0)
new_row8 <- c(8, 8, 0) 
new_row9 <- c(9, 9, 0) 
new_row10 <- c(10, 10, 0)
new_row11 <- c(11, 11, 0)
SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10, new_row11)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","grain.size","salinity"))
pred <- names(data[,env])
d <- data

################################# SERr

dcomm <- SERr2

# Dissimilarity matrix
dcomm <- cbind(year = paste(d$year), as.data.frame(dcomm))  # Adding site column for gdm function

dim(d)
dim(dcomm)


### Model ----------------------------------------------------------------------------------

# Data gdm format
print("formatting gdm data... ")
print(Sys.time())
gdm_data <- formatsitepair(bioData = dcomm, bioFormat=3, 
                           # bioData = d[,c(1,sp)], bioFormat=1, dist="bray", abundance=TRUE,
                           siteColumn="year", XColumn="long", YColumn="lat", 
                           predData = d[,c("year","long","lat", pred)])

gdm.1 <- gdm(data = gdm_data, geo = FALSE, splines = NULL, knots = NULL)

summary(gdm.1)
plot(gdm.1, plot.layout=c(2,3))

################################# SERa

dcomm <- SERa2

# Dissimilarity matrix
dcomm <- cbind(year = paste(d$year), as.data.frame(dcomm))  # Adding site column for gdm function

dim(d)
dim(dcomm)


### Model ----------------------------------------------------------------------------------

# Data gdm format
print("formatting gdm data... ")
print(Sys.time())
gdm_data <- formatsitepair(bioData = dcomm, bioFormat=3, 
                           # bioData = d[,c(1,sp)], bioFormat=1, dist="bray", abundance=TRUE,
                           siteColumn="year", XColumn="long", YColumn="lat", 
                           predData = d[,c("year","long","lat", pred)])

gdm.1 <- gdm(data = gdm_data, geo = FALSE, splines = NULL, knots = NULL)

summary(gdm.1)
plot(gdm.1, plot.layout=c(2,3))

#####################################

data <- read.csv("A_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 7:N;


#richness based turnover index  
SERr = turnover(data[, SpecColumns], method = "SERr") #explicit
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)
new_row5 <- c(5, 5, 0)
new_row6 <- c(6, 6, 0)
new_row7 <- c(7, 7, 0)
new_row8 <- c(8, 8, 0) 
new_row9 <- c(9, 9, 0) 
new_row10 <- c(10, 10, 0)
new_row11 <- c(11, 11, 0)
SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10, new_row11)
SERr2 <- acast(SERr, To ~ From , value.var="SER", na.rm= TRUE)


#abundance based turnover index 
SERa = turnover_s(data[, SpecColumns], method = "SERa")
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)
new_row5 <- c(5, 5, 0)
new_row6 <- c(6, 6, 0)
new_row7 <- c(7, 7, 0)
new_row8 <- c(8, 8, 0) 
new_row9 <- c(9, 9, 0) 
new_row10 <- c(10, 10, 0)
new_row11 <- c(11, 11, 0)
SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10, new_row11)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","grain.size","salinity"))
pred <- names(data[,env])
d <- data

################################# SERr

dcomm <- SERr2

# Dissimilarity matrix
dcomm <- cbind(year = paste(d$year), as.data.frame(dcomm))  # Adding site column for gdm function

dim(d)
dim(dcomm)


### Model ----------------------------------------------------------------------------------

# Data gdm format
print("formatting gdm data... ")
print(Sys.time())
gdm_data <- formatsitepair(bioData = dcomm, bioFormat=3, 
                           # bioData = d[,c(1,sp)], bioFormat=1, dist="bray", abundance=TRUE,
                           siteColumn="year", XColumn="long", YColumn="lat", 
                           predData = d[,c("year","long","lat", pred)])

gdm.1 <- gdm(data = gdm_data, geo = FALSE, splines = NULL, knots = NULL)

summary(gdm.1)
plot(gdm.1, plot.layout=c(2,3))

################################# SERa

dcomm <- SERa2

# Dissimilarity matrix
dcomm <- cbind(year = paste(d$year), as.data.frame(dcomm))  # Adding site column for gdm function

dim(d)
dim(dcomm)


### Model ----------------------------------------------------------------------------------

# Data gdm format
print("formatting gdm data... ")
print(Sys.time())
gdm_data <- formatsitepair(bioData = dcomm, bioFormat=3, 
                           # bioData = d[,c(1,sp)], bioFormat=1, dist="bray", abundance=TRUE,
                           siteColumn="year", XColumn="long", YColumn="lat", 
                           predData = d[,c("year","long","lat", pred)])

gdm.1 <- gdm(data = gdm_data, geo = FALSE, splines = NULL, knots = NULL)

summary(gdm.1)
plot(gdm.1, plot.layout=c(2,3))

#####################################

data <- read.csv("B_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 7:N;


#richness based turnover index  
SERr = turnover(data[, SpecColumns], method = "SERr") #explicit
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)
new_row5 <- c(5, 5, 0)
new_row6 <- c(6, 6, 0)
new_row7 <- c(7, 7, 0)
new_row8 <- c(8, 8, 0) 
new_row9 <- c(9, 9, 0) 
new_row10 <- c(10, 10, 0)
new_row11 <- c(11, 11, 0)
SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10, new_row11)
SERr2 <- acast(SERr, To ~ From , value.var="SER", na.rm= TRUE)


#abundance based turnover index 
SERa = turnover_s(data[, SpecColumns], method = "SERa")
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)
new_row5 <- c(5, 5, 0)
new_row6 <- c(6, 6, 0)
new_row7 <- c(7, 7, 0)
new_row8 <- c(8, 8, 0) 
new_row9 <- c(9, 9, 0) 
new_row10 <- c(10, 10, 0)
new_row11 <- c(11, 11, 0)
SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10, new_row11)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","grain.size","salinity"))
pred <- names(data[,env])
d <- data

################################# SERr

dcomm <- SERr2

# Dissimilarity matrix
dcomm <- cbind(year = paste(d$year), as.data.frame(dcomm))  # Adding site column for gdm function

dim(d)
dim(dcomm)


### Model ----------------------------------------------------------------------------------

# Data gdm format
print("formatting gdm data... ")
print(Sys.time())
gdm_data <- formatsitepair(bioData = dcomm, bioFormat=3, 
                           # bioData = d[,c(1,sp)], bioFormat=1, dist="bray", abundance=TRUE,
                           siteColumn="year", XColumn="long", YColumn="lat", 
                           predData = d[,c("year","long","lat", pred)])

gdm.1 <- gdm(data = gdm_data, geo = FALSE, splines = NULL, knots = NULL)

summary(gdm.1)
plot(gdm.1, plot.layout=c(2,3))

################################# SERa

dcomm <- SERa2

# Dissimilarity matrix
dcomm <- cbind(year = paste(d$year), as.data.frame(dcomm))  # Adding site column for gdm function

dim(d)
dim(dcomm)


### Model ----------------------------------------------------------------------------------

# Data gdm format
print("formatting gdm data... ")
print(Sys.time())
gdm_data <- formatsitepair(bioData = dcomm, bioFormat=3, 
                           # bioData = d[,c(1,sp)], bioFormat=1, dist="bray", abundance=TRUE,
                           siteColumn="year", XColumn="long", YColumn="lat", 
                           predData = d[,c("year","long","lat", pred)])

gdm.1 <- gdm(data = gdm_data, geo = FALSE, splines = NULL, knots = NULL)

summary(gdm.1)
plot(gdm.1, plot.layout=c(2,3))

#####################################

data <- read.csv("M_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 7:N;


#richness based turnover index  
SERr = turnover(data[, SpecColumns], method = "SERr") #explicit
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)
new_row5 <- c(5, 5, 0)
new_row6 <- c(6, 6, 0)
new_row7 <- c(7, 7, 0)
new_row8 <- c(8, 8, 0) 
new_row9 <- c(9, 9, 0) 
new_row10 <- c(10, 10, 0)
new_row11 <- c(11, 11, 0)
SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10, new_row11)
SERr2 <- acast(SERr, To ~ From , value.var="SER", na.rm= TRUE)


#abundance based turnover index 
SERa = turnover_s(data[, SpecColumns], method = "SERa")
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)
new_row5 <- c(5, 5, 0)
new_row6 <- c(6, 6, 0)
new_row7 <- c(7, 7, 0)
new_row8 <- c(8, 8, 0) 
new_row9 <- c(9, 9, 0) 
new_row10 <- c(10, 10, 0)
new_row11 <- c(11, 11, 0)
SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10, new_row11)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","grain.size","salinity"))
pred <- names(data[,env])
d <- data

################################# SERr

dcomm <- SERr2

# Dissimilarity matrix
dcomm <- cbind(year = paste(d$year), as.data.frame(dcomm))  # Adding site column for gdm function

dim(d)
dim(dcomm)


### Model ----------------------------------------------------------------------------------

# Data gdm format
print("formatting gdm data... ")
print(Sys.time())
gdm_data <- formatsitepair(bioData = dcomm, bioFormat=3, 
                           # bioData = d[,c(1,sp)], bioFormat=1, dist="bray", abundance=TRUE,
                           siteColumn="year", XColumn="long", YColumn="lat", 
                           predData = d[,c("year","long","lat", pred)])

gdm.1 <- gdm(data = gdm_data, geo = FALSE, splines = NULL, knots = NULL)

summary(gdm.1)
plot(gdm.1, plot.layout=c(2,3))

################################# SERa

dcomm <- SERa2

# Dissimilarity matrix
dcomm <- cbind(year = paste(d$year), as.data.frame(dcomm))  # Adding site column for gdm function

dim(d)
dim(dcomm)


### Model ----------------------------------------------------------------------------------

# Data gdm format
print("formatting gdm data... ")
print(Sys.time())
gdm_data <- formatsitepair(bioData = dcomm, bioFormat=3, 
                           # bioData = d[,c(1,sp)], bioFormat=1, dist="bray", abundance=TRUE,
                           siteColumn="year", XColumn="long", YColumn="lat", 
                           predData = d[,c("year","long","lat", pred)])

gdm.1 <- gdm(data = gdm_data, geo = FALSE, splines = NULL, knots = NULL)

summary(gdm.1)
plot(gdm.1, plot.layout=c(2,3))

