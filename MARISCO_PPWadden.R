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

ppwadden <- read.csv("MARISCO_pp_wadden_yearly_StandingDiv_Abu.csv")

ppwadden$StationID <- as.factor(ppwadden$StationID)

#Are the missing values?

levels(ppwadden$StationID)
colSums(is.na(ppwadden))

ppwadden <- na.exclude(ppwadden) 

dim(ppwadden)

levels(ppwadden$StationID)

summary(ppwadden)

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


boxplot(ppwadden$biomass, 
        main = "Carbon-based biomass")
hist(ppwadden$biomass)
dotchart(ppwadden$biomass)

par(mfrow = c(1, 2))
plot(x = ppwadden$year, 
     y = ppwadden$S)
plot(x = ppwadden$year, 
     y = ppwadden$Sy)
plot(x = ppwadden$year, 
     y = ppwadden$ENS)
plot(x = ppwadden$year, 
     y = ppwadden$ENSy)
par(mfrow = c(1, 3))
plot(x = ppwadden$year, 
     y = ppwadden$biomass)
par(mfrow = c(1, 1))

# Collinearity X

names(ppwadden)
pairs.panels(ppwadden[, c("salinity", "temperature", "total.n", "total.p", "pH", 
                          "S", "ENS", "Sy", "ENSy", "biomass")])

# Zero inflation
sum(ppwadden == 0)

# Scale variables (method for GLMs)
ppwadden$salinity <- scale(ppwadden$salinity)
ppwadden$temperature <- scale(ppwadden$temperature)
ppwadden$total.n <- scale(ppwadden$total.n)
ppwadden$pH <- scale(ppwadden$pH)
ppwadden$total.p <- scale(ppwadden$total.p)

# Recode vars to roughly same scale
#ppwadden$biomass <- ppwadden$mean.bm/10000
ppwadden$biomass <- log(ppwadden$biomass)
ppwadden$S <- ppwadden$S/10
ppwadden$Sy <- ppwadden$Sy/10



# ################################# gls fitting


gls.abu.S <- gls(S ~ temperature + salinity + total.p + total.n, data=ppwadden,
             correlation=corGaus(form= ~ lat + long | year))
summary(gls.abu.S)
plot(gls.abu.S)
hist(resid(gls.abu.S,type="pearson"))

gls.abu.ENS <- gls(ENS ~ temperature + salinity + total.p + total.n, data=ppwadden,
               correlation=corGaus(form= ~ long + lat | year))
summary(gls.abu.ENS)
plot(gls.abu.ENS)
hist(resid(gls.abu.ENS,type="pearson"))

gls.abu.biom <- gls(biomass ~ temperature + salinity + pH + total.n + total.p + ENS + S, data=ppwadden,
                correlation=corGaus(form= ~ long + lat | year))
summary(gls.abu.biom)
plot(gls.abu.biom)
hist(resid(gls.abu.biom,type="pearson"))

gls.abu.Sy <- gls(Sy ~ temperature + salinity + pH + total.n, data=ppwadden,
                  correlation=corGaus(form= ~ lat + long | year))
summary(gls.abu.Sy)
plot(gls.abu.Sy)
hist(resid(gls.abu.Sy,type="pearson"))

gls.abu.ENSy <- gls(ENSy ~ temperature + salinity + pH + total.p, data=ppwadden,
                    correlation=corGaus(form= ~ long + lat | year))
summary(gls.abu.ENSy)
plot(gls.abu.ENSy)
hist(resid(gls.abu.ENSy,type="pearson"))

gls.abu.biomy <- gls(biomass ~ temperature + salinity + pH + total.n + total.p + ENSy + Sy, data=ppwadden,
                 correlation=corGaus(form= ~ long + lat | year))
summary(gls.abu.biomy)
plot(gls.abu.biomy)
hist(resid(gls.abu.biomy,type="pearson"))


################################# piecewiseSEM fitting


sem.SD.abu.gls1 <- psem(
                gls.abu.S,
                gls.abu.ENS,
                gls.abu.biom,
                S %~~% ENS,
                total.n %~~% total.p,
                total.n %~~% salinity,
                total.p %~~% salinity,
                pH %~~% biomass,
                ppwadden)
summary(sem.SD.abu.gls1)

sem.SD.abu.gls2 <- psem(
                        gls.abu.Sy,
                        gls.abu.ENSy,
                        gls.abu.biomy,
                        Sy %~~% ENSy,
                        total.n %~~% total.p,
                        total.n %~~% salinity,
                        total.p %~~% salinity,
                        pH %~~% biomass,
                        ppwadden)
summary(sem.SD.abu.gls2)
summary(sem.SD.abu.gls1)
plot(sem.SD.abu.gls1)
plot(sem.SD.abu.gls2)

################################################## GDM Model

source("turnover.R")

#####################################

data <- read.csv("ZUIDOLWOT_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 9:N;


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
SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8)
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
SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","pH","salinity", "total.p","total.n"))
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

data <- read.csv("MARSDND_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 9:N;


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
SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8)
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
SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","pH","salinity", "total.p","total.n"))
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
SpecColumns = 9:N;


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
SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8)
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
SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","pH","salinity", "total.p","total.n"))
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

data <- read.csv("TERSLG4_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 9:N;


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
env <- which(colnames(data) %in% c("temperature","pH","salinity", "total.p","total.n"))
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

data <- read.csv("GROOTGND_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 9:N;


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
SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8)
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
SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","pH","salinity", "total.p","total.n"))
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
SpecColumns = 9:N;


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
SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8)
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
SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","pH","salinity", "total.p","total.n"))
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
SpecColumns = 9:N;


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
SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8)
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
SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","pH","salinity", "total.p","total.n"))
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
SpecColumns = 9:N;


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
env <- which(colnames(data) %in% c("temperature","pH","salinity", "total.p","total.n"))
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

# Scale variables (method for GLMs)
data$salinity <- scale(data$salinity)
data$temperature <- scale(data$temperature)
data$total.n <- scale(data$total.n)
data$pH <- scale(data$pH)
data$total.p <- scale(data$total.p)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 9:N;


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
SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","pH","salinity", "total.p","total.n"))
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
SpecColumns = 9:N;

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
env <- which(colnames(data) %in% c("temperature","pH","salinity", "total.p","total.n"))
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
SpecColumns = 9:N;

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
SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","pH","salinity", "total.p","total.n"))
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

data <- read.csv("selected_pp_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

# Scale variables (method for GLMs)
data$salinity <- scale(data$salinity)
data$temperature <- scale(data$temperature)
data$total.n <- scale(data$total.n)
data$pH <- scale(data$pH)
data$total.p <- scale(data$total.p)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 9:N;

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
SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, 
              new_row10, new_row11, new_row12, new_row13, new_row14, new_row15, new_row16, new_row17, new_row18, 
              new_row19)
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
new_row13 <- c(13, 13, 0)
new_row14 <- c(14, 14, 0)
new_row15 <- c(15, 15, 0)
new_row16 <- c(16, 16, 0)
new_row17 <- c(17, 17, 0)
new_row18 <- c(18, 18, 0)
new_row19 <- c(19, 19, 0)
SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, 
              new_row10, new_row11, new_row12, new_row13, new_row14, new_row15, new_row16, new_row17, new_row18, 
              new_row19)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","pH","salinity", "total.p","total.n"))
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

