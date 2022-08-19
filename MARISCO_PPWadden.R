setwd(".../Data/Phytoplankton")

######################################################################

library(lattice)
source("HighstatLibV10.R")
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
library(reshape2)
library(gdm)

###################################################################### Standing Diversity Model
####################################################### SEM All stations

ppwadden <- read.csv("INTERREG_pp_wadden_annual.csv")

ppwadden$StationID <- as.factor(ppwadden$StationID)

#Are the missing values?

levels(ppwadden$StationID)
colSums(is.na(ppwadden))

ppwadden <- na.exclude(ppwadden) 

dim(ppwadden)

# Outliers

par(mfrow = c(1, 3))
boxplot(ppwadden$S, 
        main = "Species richness abundance")
hist(ppwadden$S)
dotchart(ppwadden$S)
boxplot(ppwadden$ENS, 
        main = "Effective number of species abundance")
hist(ppwadden$ENS)
dotchart(ppwadden$ENS)
boxplot(ppwadden$Sy, 
        main = "Species richness y abundance")
hist(ppwadden$Sy)
dotchart(ppwadden$Sy)
boxplot(ppwadden$ENSy, 
        main = "Effective number of species y abundance")
hist(ppwadden$ENSy)
dotchart(ppwadden$ENSy)


boxplot(ppwadden$carbon, 
        main = "Carbon-based biomass")
hist(ppwadden$carbon)
dotchart(ppwadden$carbon)

par(mfrow = c(1, 2))
plot(x = ppwadden$year, 
     y = ppwadden$S)
plot(x = ppwadden$year, 
     y = ppwadden$Sy)
plot(x = ppwadden$year, 
     y = ppwadden$ENS)
plot(x = ppwadden$year, 
     y = ppwadden$ENSy)
plot(x = ppwadden$year, 
     y = ppwadden$carbon)
par(mfrow = c(1, 1))

# Collinearity X

names(ppwadden)
pairs.panels(ppwadden[, c("salinity", "temperature", "total.n", "din", "dip", "total.p", "pH", "silicon", 
                          "suspended.particulates", "S", "ENS", "carbon")])

# Zero inflation
sum(ppwadden == 0)

# Scale variables (method for GLMs)
ppwadden$salinity <- scale(ppwadden$salinity)
ppwadden$temperature <- scale(ppwadden$temperature)
ppwadden$total.n <- scale(ppwadden$total.n)
ppwadden$total.p <- scale(ppwadden$total.p)
ppwadden$din <- scale(ppwadden$din)
ppwadden$dip <- scale(ppwadden$dip)
ppwadden$pH <- scale(ppwadden$pH)
ppwadden$silicon <- scale(ppwadden$silicon)
ppwadden$suspended.particulates <- scale(ppwadden$suspended.particulates)

# Recode vars to roughly same scale
#ppwadden$biomass <- ppwadden$mean.bm/10000
ppwadden$carbon <- log(ppwadden$carbon)
ppwadden$S <- ppwadden$S/10
ppwadden$Sy <- ppwadden$Sy/10
ppwadden$ENS <- log(ppwadden$ENS)
ppwadden$ENSy <- log(ppwadden$ENSy)

################################### gls fitting


gls.abu.S <- gls(S ~ temperature + salinity + total.p  + silicon, 
                 data=ppwadden, correlation=corGaus(form= ~ lat + long | year))
summary(gls.abu.S)
plot(gls.abu.S)
hist(resid(gls.abu.S,type="pearson"))

gls.abu.ENS <- gls(ENS ~ temperature + salinity + total.n  + silicon, 
                   data=ppwadden, correlation=corGaus(form= ~ long + lat | year))
summary(gls.abu.ENS)
plot(gls.abu.ENS)
hist(resid(gls.abu.ENS,type="pearson"))

gls.abu.biom <- gls(carbon ~ temperature + salinity + pH + total.p + total.n + silicon + ENS + S, data=ppwadden,
                    correlation=corGaus(form= ~ long + lat | year))
summary(gls.abu.biom)
plot(gls.abu.biom)
hist(resid(gls.abu.biom,type="pearson"))



################################# piecewiseSEM fitting

sem.gls.carbon <- psem(
  gls.abu.S,
  gls.abu.ENS,
  gls.abu.biom,
  S %~~% ENS,
  total.n %~~% total.p,
  total.n %~~% salinity,
  total.p %~~% salinity,
  pH %~~% carbon,
  silicon %~~% total.p,
  silicon %~~% total.n,
  silicon %~~% salinity,
  ppwadden)
summary(sem.gls.carbon)
plot(sem.gls.carbon)


################################################## GDM Model

source("turnover.R")

#####################################

data <- read.csv("MARSDND_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 10:N;


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
new_row13 <- c(13, 13, 0)
new_row14 <- c(14, 14, 0)
new_row15 <- c(15, 15, 0)
new_row16 <- c(16, 16, 0)
new_row17 <- c(17, 17, 0)
new_row18 <- c(18, 18, 0) 
new_row19 <- c(19, 19, 0)

SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10
              , new_row11, new_row12, new_row13, new_row14, new_row15, new_row16, new_row17, new_row18, new_row19)
SERr2 <- acast(SERr, To ~ From , value.var="SER", na.rm= TRUE)


#abundance based turnover index 
SERa = turnover(data[, SpecColumns], method = "SERa") #explicit
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
new_row13 <- c(13, 13, 0)
new_row14 <- c(14, 14, 0)
new_row15 <- c(15, 15, 0)
new_row16 <- c(16, 16, 0)
new_row17 <- c(17, 17, 0)
new_row18 <- c(18, 18, 0) 
new_row19 <- c(19, 19, 0)

SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10
              , new_row11, new_row12, new_row13, new_row14, new_row15, new_row16, new_row17, new_row18, new_row19)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","pH","salinity", "total.p","total.n", "silicon"))
pred <- names(data[,env])
d <- data

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

######### SERr

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
#####################################

data <- read.csv("HUIBGOT_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 10:N;


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
new_row13 <- c(13, 13, 0)
new_row14 <- c(14, 14, 0)
new_row15 <- c(15, 15, 0)
new_row16 <- c(16, 16, 0)
new_row17 <- c(17, 17, 0)
new_row18 <- c(18, 18, 0) 
new_row19 <- c(19, 19, 0)

SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10
              , new_row11, new_row12, new_row13, new_row14, new_row15, new_row16, new_row17, new_row18, new_row19)
SERr2 <- acast(SERr, To ~ From , value.var="SER", na.rm= TRUE)


#abundance based turnover index 
SERa = turnover(data[, SpecColumns], method = "SERa") #explicit
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
new_row13 <- c(13, 13, 0)
new_row14 <- c(14, 14, 0)
new_row15 <- c(15, 15, 0)
new_row16 <- c(16, 16, 0)
new_row17 <- c(17, 17, 0)
new_row18 <- c(18, 18, 0) 
new_row19 <- c(19, 19, 0)

SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10
              , new_row11, new_row12, new_row13, new_row14, new_row15, new_row16, new_row17, new_row18, new_row19)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","pH","salinity", "total.p","total.n", "silicon"))
pred <- names(data[,env])
d <- data

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

######### SERr

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

###############################################

data <- read.csv("GROOTGND_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 10:N;


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
new_row13 <- c(13, 13, 0)
new_row14 <- c(14, 14, 0)
new_row15 <- c(15, 15, 0)
new_row16 <- c(16, 16, 0)
new_row17 <- c(17, 17, 0)
new_row18 <- c(18, 18, 0) 
new_row19 <- c(19, 19, 0)

SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10
              , new_row11, new_row12, new_row13, new_row14, new_row15, new_row16, new_row17, new_row18, new_row19)
SERr2 <- acast(SERr, To ~ From , value.var="SER", na.rm= TRUE)


#abundance based turnover index 
SERa = turnover(data[, SpecColumns], method = "SERa") #explicit
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
new_row13 <- c(13, 13, 0)
new_row14 <- c(14, 14, 0)
new_row15 <- c(15, 15, 0)
new_row16 <- c(16, 16, 0)
new_row17 <- c(17, 17, 0)
new_row18 <- c(18, 18, 0) 
new_row19 <- c(19, 19, 0)

SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10
              , new_row11, new_row12, new_row13, new_row14, new_row15, new_row16, new_row17, new_row18, new_row19)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","pH","salinity", "total.p","total.n", "silicon"))
pred <- names(data[,env])
d <- data

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

######### SERr

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


#####################################

data <- read.csv("DANTZGT_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 10:N;


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
new_row13 <- c(13, 13, 0)
new_row14 <- c(14, 14, 0)
new_row15 <- c(15, 15, 0)
new_row16 <- c(16, 16, 0)
new_row17 <- c(17, 17, 0)
new_row18 <- c(18, 18, 0) 
new_row19 <- c(19, 19, 0)

SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10
              , new_row11, new_row12, new_row13, new_row14, new_row15, new_row16, new_row17, new_row18, new_row19)
SERr2 <- acast(SERr, To ~ From , value.var="SER", na.rm= TRUE)


#abundance based turnover index 
SERa = turnover(data[, SpecColumns], method = "SERa") #explicit
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
new_row13 <- c(13, 13, 0)
new_row14 <- c(14, 14, 0)
new_row15 <- c(15, 15, 0)
new_row16 <- c(16, 16, 0)
new_row17 <- c(17, 17, 0)
new_row18 <- c(18, 18, 0) 
new_row19 <- c(19, 19, 0)

SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10
              , new_row11, new_row12, new_row13, new_row14, new_row15, new_row16, new_row17, new_row18, new_row19)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","pH","salinity", "total.p","total.n", "silicon"))
pred <- names(data[,env])
d <- data

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

######### SERr

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


#####################################

data <- read.csv("ROTTMPT3_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 10:N;


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
new_row13 <- c(13, 13, 0)
new_row14 <- c(14, 14, 0)
new_row15 <- c(15, 15, 0)
new_row16 <- c(16, 16, 0)
new_row17 <- c(17, 17, 0)
new_row18 <- c(18, 18, 0) 
new_row19 <- c(19, 19, 0)

SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10
              , new_row11, new_row12, new_row13, new_row14, new_row15, new_row16, new_row17, new_row18, new_row19)
SERr2 <- acast(SERr, To ~ From , value.var="SER", na.rm= TRUE)


#abundance based turnover index 
SERa = turnover(data[, SpecColumns], method = "SERa") #explicit
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
new_row13 <- c(13, 13, 0)
new_row14 <- c(14, 14, 0)
new_row15 <- c(15, 15, 0)
new_row16 <- c(16, 16, 0)
new_row17 <- c(17, 17, 0)
new_row18 <- c(18, 18, 0) 
new_row19 <- c(19, 19, 0)

SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10
              , new_row11, new_row12, new_row13, new_row14, new_row15, new_row16, new_row17, new_row18, new_row19)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","pH","salinity", "total.p","total.n", "silicon"))
pred <- names(data[,env])
d <- data

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

######### SERr

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


#####################################

data <- read.csv("BOOMKDP_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 10:N;


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

SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10
              , new_row11, new_row12)
SERr2 <- acast(SERr, To ~ From , value.var="SER", na.rm= TRUE)


#abundance based turnover index 
SERa = turnover(data[, SpecColumns], method = "SERa") #explicit
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

SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10
              , new_row11, new_row12)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","pH","salinity", "total.p","total.n", "silicon"))
pred <- names(data[,env])
d <- data

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

######### SERr

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

#####################################

data <- read.csv("DOOVBWT_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 10:N;


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

SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10
              , new_row11, new_row12)
SERr2 <- acast(SERr, To ~ From , value.var="SER", na.rm= TRUE)


#abundance based turnover index 
SERa = turnover(data[, SpecColumns], method = "SERa") #explicit
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

SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10
              , new_row11, new_row12)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","pH","salinity", "total.p","total.n", "silicon"))
pred <- names(data[,env])
d <- data

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

######### SERr

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

#####################################

data <- read.csv("TERSLG10_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 10:N;


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

SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10
              , new_row11, new_row12)
SERr2 <- acast(SERr, To ~ From , value.var="SER", na.rm= TRUE)


#abundance based turnover index 
SERa = turnover(data[, SpecColumns], method = "SERa") #explicit
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

SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10
              , new_row11, new_row12)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","pH","salinity", "total.p","total.n", "silicon"))
pred <- names(data[,env])
d <- data

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

######### SERr

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


#####################################

data <- read.csv("Bork_W_1_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 10:N;


#richness based turnover index  
SERr = turnover(data[, SpecColumns], method = "SERr") #explicit
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)
new_row5 <- c(5, 5, 0)
new_row6 <- c(6, 6, 0)

SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6)
SERr2 <- acast(SERr, To ~ From , value.var="SER", na.rm= TRUE)


#abundance based turnover index 
SERa = turnover(data[, SpecColumns], method = "SERa") #explicit
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)
new_row5 <- c(5, 5, 0)
new_row6 <- c(6, 6, 0)

SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","pH","salinity", "total.p","total.n", "silicon"))
pred <- names(data[,env])
d <- data

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

######### SERr

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


#####################################

data <- read.csv("JaBu_W_1_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 10:N;


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

SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10)
SERr2 <- acast(SERr, To ~ From , value.var="SER", na.rm= TRUE)


#abundance based turnover index 
SERa = turnover(data[, SpecColumns], method = "SERa") #explicit
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

SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","pH","salinity", "total.p","total.n", "silicon"))
pred <- names(data[,env])
d <- data

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

######### SERr

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


#####################################

data <- read.csv("Nney_W_2_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 10:N;


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

SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9,
              new_row10, new_row11)
SERr2 <- acast(SERr, To ~ From , value.var="SER", na.rm= TRUE)


#abundance based turnover index 
SERa = turnover(data[, SpecColumns], method = "SERa") #explicit
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

SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9,
              new_row10, new_row11)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","total.p","total.n", "silicon"))
pred <- names(data[,env])
d <- data

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

######### SERr

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

#####################################

data <- read.csv("WeMu_W_1_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 10:N;


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
SERa = turnover(data[, SpecColumns], method = "SERa") #explicit
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
env <- which(colnames(data) %in% c("temperature","pH","salinity", "total.p","total.n", "silicon"))
pred <- names(data[,env])
d <- data

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

######### SERr

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

