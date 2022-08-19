setwd(".../Birds")

library(lattice)
library(dplyr)
library(mgcv)
library(piecewiseSEM)
library(lme4)
library(psych)
library(plyr)
library(tidyverse)
library(reshape2)
library(lubridate)
library(vegan)
library(sjPlot)
library(MASS)
library(nlme)
library(sjPlot)
library(ggplot2)
library(gdm)
library(tidyverse)
library(vegan)

####################################################### SEM Analyis

data <- read.csv("MARISCO_birds_wadden_SEM.csv")


data$StationID <- as.factor(data$StationID)
unique(data$year)
unique(data$StationID)

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

summary(data)

pairs.panels(data[,2:11])


# Outliers

par(mfrow = c(1, 3))
boxplot(data$S, 
        main = "Species richness abundance")
hist(data$S)
dotchart(data$S)
boxplot(data$ENS, 
        main = "Effective number of species abundance")
hist(data$ENS)
dotchart(data$ENS)
boxplot(data$Sy, 
        main = "Species richness y abundance")
hist(data$Sy)
dotchart(data$Sy)
boxplot(data$ENSy, 
        main = "Effective number of species y abundance")
hist(data$ENSy)
dotchart(data$ENSy)

par(mfrow = c(1, 2))
plot(x = data$year, 
     y = data$S)
plot(x = data$year, 
     y = data$Sy)
plot(x = data$year, 
     y = data$ENS)
plot(x = data$year, 
     y = data$ENSy)

par(mfrow = c(1, 1))

# Zero inflation
sum(data == 0)

# Scale variables (method for GLMs)
data$precipitation <- scale(data$precipitation)
data$temperature <- scale(data$temperature)
data$mzb.biomass <- scale(data$mzb.biomass)

# Recode vars to roughly same scale

data$S <- data$S/10
data$Sy <- data$Sy/10

################################# gls fitting 


gls.abu.S <- gls(S ~ temperature + precipitation + mzb.biomass, data=data, 
                 correlation=corGaus(form= ~ lat + long | year))
summary(gls.abu.S)
plot(gls.abu.S)

gls.abu.ENS <- gls(ENS ~  temperature + precipitation + mzb.biomass, data=data, 
                   correlation=corGaus(form= ~ long + lat | year))
summary(gls.abu.ENS)
plot(gls.abu.ENS)

gls.abu.Sy <- gls(Sy ~  temperature + precipitation + mzb.biomass, data=data, 
                  correlation=corGaus(form= ~ lat + long | year))
summary(gls.abu.Sy)
plot(gls.abu.Sy)

gls.abu.ENSy <- gls(ENSy ~  temperature + precipitation + mzb.biomass, data=data, 
                    correlation=corGaus(form= ~ long + lat | year))
summary(gls.abu.ENSy)
plot(gls.abu.ENSy)


################################# piecewiseSEM fitting

sem.SD.abu.gls1 <- psem(gls.abu.S,
                        gls.abu.ENS,
                        S %~~% ENS,
                        data)
summary(sem.SD.abu.gls1)

sem.SD.abu.gls2 <- psem(gls.abu.Sy,
                        gls.abu.ENSy,
                        Sy %~~% ENSy,
                        data)
summary(sem.SD.abu.gls2)
summary(sem.SD.abu.gls1)


################################################## GDM Model

source("turnover.R")

data <- read.csv("leybucht_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)


# Scale variables (method for GLMs)
data$precipitation <- scale(data$precipitation)
data$temperature <- scale(data$temperature)
data$mzb.biomass <- scale(data$mzb.biomass)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 7:N;


# #richness based turnover index  
SERr = turnover_s(data[, SpecColumns], method = "SERr") #explicit
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
SERa = turnover_s(data[, SpecColumns], method = "SERa") #explicit
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



################################ Code Marina

sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","precipitation","mzb.biomass"))
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


################ SERr


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


#########################################################################

data <- read.csv("rysumer_nacken_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)


# Scale variables (method for GLMs)
data$precipitation <- scale(data$precipitation)
data$temperature <- scale(data$temperature)
data$mzb.biomass <- scale(data$mzb.biomass)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 7:N;


# #richness based turnover index  
SERr = turnover_s(data[, SpecColumns], method = "SERr") #explicit
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)

SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4)
SERr2 <- acast(SERr, To ~ From , value.var="SER", na.rm= TRUE)

#abundance based turnover index 
SERa = turnover_s(data[, SpecColumns], method = "SERa") #explicit
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)

SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)



################################ Code Marina

sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","precipitation","mzb.biomass"))
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


################ SERr


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


#########################################################################

data <- read.csv("borkum_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)


# Scale variables (method for GLMs)
data$precipitation <- scale(data$precipitation)
data$temperature <- scale(data$temperature)
data$mzb.biomass <- scale(data$mzb.biomass)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 7:N;

# #richness based turnover index  
SERr = turnover_s(data[, SpecColumns], method = "SERr") #explicit
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
SERa = turnover_s(data[, SpecColumns], method = "SERa") #explicit
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

################################ Code Marina

sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","precipitation","mzb.biomass"))
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


################ SERr


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



#########################################################################

data <- read.csv("hilgenriedersiel_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)


# Scale variables (method for GLMs)
data$precipitation <- scale(data$precipitation)
data$temperature <- scale(data$temperature)
data$mzb.biomass <- scale(data$mzb.biomass)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 7:N;

# #richness based turnover index  
SERr = turnover_s(data[, SpecColumns], method = "SERr") #explicit
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

SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9,
              new_row10, new_row11, new_row12, new_row13, new_row14, new_row15, new_row16, new_row17, new_row18)
SERr2 <- acast(SERr, To ~ From , value.var="SER", na.rm= TRUE)

#abundance based turnover index 
SERa = turnover_s(data[, SpecColumns], method = "SERa") #explicit
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

SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, 
              new_row10, new_row11, new_row12, new_row13, new_row14, new_row15, new_row16, new_row17, new_row18)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)

################################ Code Marina

sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","precipitation","mzb.biomass"))
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


################ SERr


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

###############################################################

data <- read.csv("norderney_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

# Scale variables (method for GLMs)
data$precipitation <- scale(data$precipitation)
data$temperature <- scale(data$temperature)
data$mzb.biomass <- scale(data$mzb.biomass)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 7:N;


# #richness based turnover index  
SERr = turnover_s(data[, SpecColumns], method = "SERr") #explicit
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
new_row20 <- c(20, 20, 0)
new_row21 <- c(21, 21, 0)
new_row22 <- c(22, 22, 0)
new_row23 <- c(23, 23, 0)
new_row24 <- c(24, 24, 0)
new_row25 <- c(25, 25, 0)
new_row26 <- c(26, 26, 0)
new_row27 <- c(27, 27, 0)
new_row28 <- c(28, 28, 0)
new_row29 <- c(29, 29, 0)
new_row30 <- c(30, 30, 0)
new_row31 <- c(31, 31, 0)
new_row32 <- c(32, 32, 0)
new_row33 <- c(33, 33, 0)
new_row34 <- c(34, 34, 0)
new_row35 <- c(35, 35, 0)
new_row36 <- c(36, 36, 0)
SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9,
              new_row10, new_row11, new_row12, new_row13, new_row14, new_row15, new_row16, new_row17, new_row18,
              new_row19, new_row20, new_row21, new_row22, new_row23, new_row24, new_row25, new_row26, new_row27,
              new_row28, new_row29, new_row30, new_row31, new_row32, new_row33, new_row34, new_row35, new_row36)
SERr2 <- acast(SERr, To ~ From , value.var="SER", na.rm= TRUE)

#abundance based turnover index 
SERa = turnover_s(data[, SpecColumns], method = "SERa") #explicit
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
new_row20 <- c(20, 20, 0)
new_row21 <- c(21, 21, 0)
new_row22 <- c(22, 22, 0)
new_row23 <- c(23, 23, 0)
new_row24 <- c(24, 24, 0)
new_row25 <- c(25, 25, 0)
new_row26 <- c(26, 26, 0)
new_row27 <- c(27, 27, 0)
new_row28 <- c(28, 28, 0)
new_row29 <- c(29, 29, 0)
new_row30 <- c(30, 30, 0)
new_row31 <- c(31, 31, 0)
new_row32 <- c(32, 32, 0)
new_row33 <- c(33, 33, 0)
new_row34 <- c(34, 34, 0)
new_row35 <- c(35, 35, 0)
new_row36 <- c(36, 36, 0)
SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, 
              new_row10, new_row11, new_row12, new_row13, new_row14, new_row15, new_row16, new_row17, new_row18, 
              new_row19, new_row20, new_row21, new_row22, new_row23, new_row24, new_row25, new_row26, new_row27, 
              new_row28, new_row29, new_row30, new_row31, new_row32, new_row33, new_row34, new_row35, new_row36)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)



################################ Code Marina

sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","precipitation","mzb.biomass"))
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

################ SERr


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


###############################################################

data <- read.csv("pilsum_manslagt_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

# Scale variables (method for GLMs)
data$precipitation <- scale(data$precipitation)
data$temperature <- scale(data$temperature)
data$mzb.biomass <- scale(data$mzb.biomass)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 7:N;



# #richness based turnover index  
SERr = turnover_s(data[, SpecColumns], method = "SERr") #explicit
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)
new_row5 <- c(5, 5, 0)
new_row6 <- c(6, 6, 0)

SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6)
SERr2 <- acast(SERr, To ~ From , value.var="SER", na.rm= TRUE)

#abundance based turnover index 
SERa = turnover_s(data[, SpecColumns], method = "SERa") #explicit
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)
new_row5 <- c(5, 5, 0)
new_row6 <- c(6, 6, 0)

SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)



################################ Code Marina

sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","precipitation","mzb.biomass"))
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

################ SERr


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

###############################################################

data <- read.csv("bremerhaven_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

# Scale variables (method for GLMs)
data$precipitation <- scale(data$precipitation)
data$temperature <- scale(data$temperature)
data$mzb.biomass <- scale(data$mzb.biomass)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 7:N;



# #richness based turnover index  
SERr = turnover_s(data[, SpecColumns], method = "SERr") #explicit
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

SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9,
              new_row10, new_row11, new_row12, new_row13, new_row14)
SERr2 <- acast(SERr, To ~ From , value.var="SER", na.rm= TRUE)

#abundance based turnover index 
SERa = turnover_s(data[, SpecColumns], method = "SERa") #explicit
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
SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, 
              new_row10, new_row11, new_row12, new_row13, new_row14)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


################################ Code Marina

sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","precipitation","mzb.biomass"))
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

################ SERr


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


###############################################################

data <- read.csv("jadebusen_ost_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

# Scale variables (method for GLMs)
data$precipitation <- scale(data$precipitation)
data$temperature <- scale(data$temperature)
data$mzb.biomass <- scale(data$mzb.biomass)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 7:N;



# #richness based turnover index  
SERr = turnover_s(data[, SpecColumns], method = "SERr") #explicit
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
SERa = turnover_s(data[, SpecColumns], method = "SERa") #explicit
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


################################ Code Marina

sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","precipitation","mzb.biomass"))
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

################ SERr


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



###############################################################

data <- read.csv("jadebusen_west_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

# Scale variables (method for GLMs)
data$precipitation <- scale(data$precipitation)
data$temperature <- scale(data$temperature)
data$mzb.biomass <- scale(data$mzb.biomass)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 7:N;



# #richness based turnover index  
SERr = turnover_s(data[, SpecColumns], method = "SERr") #explicit
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
SERa = turnover_s(data[, SpecColumns], method = "SERa") #explicit
new_row1 <- c(1, 1, 0)
new_row2 <- c(2, 2, 0)
new_row3 <- c(3, 3, 0)
new_row4 <- c(4, 4, 0)
new_row5 <- c(5, 5, 0)
new_row6 <- c(6, 6, 0)
new_row7 <- c(7, 7, 0)

SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


################################ Code Marina

sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","precipitation","mzb.biomass"))
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

################ SERr


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



###############################################################

data <- read.csv("langenwarder_deich_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

# Scale variables (method for GLMs)
data$precipitation <- scale(data$precipitation)
data$temperature <- scale(data$temperature)
data$mzb.biomass <- scale(data$mzb.biomass)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 7:N;



# #richness based turnover index  
SERr = turnover_s(data[, SpecColumns], method = "SERr") #explicit
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
SERa = turnover_s(data[, SpecColumns], method = "SERa") #explicit
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


################################ Code Marina

sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","precipitation","mzb.biomass"))
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

################ SERr


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




###############################################################

data <- read.csv("langluetjen_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

# Scale variables (method for GLMs)
data$precipitation <- scale(data$precipitation)
data$temperature <- scale(data$temperature)
data$mzb.biomass <- scale(data$mzb.biomass)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 7:N;



# #richness based turnover index  
SERr = turnover_s(data[, SpecColumns], method = "SERr") #explicit
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
new_row20 <- c(20, 20, 0)
new_row21 <- c(21, 21, 0)

SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9,
              new_row10, new_row11, new_row12, new_row13, new_row14, new_row15, new_row16, new_row17, new_row18,
              new_row19, new_row20, new_row21)
SERr2 <- acast(SERr, To ~ From , value.var="SER", na.rm= TRUE)

#abundance based turnover index 
SERa = turnover_s(data[, SpecColumns], method = "SERa") #explicit
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
new_row20 <- c(20, 20, 0)
new_row21 <- c(21, 21, 0)

SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, 
              new_row10, new_row11, new_row12, new_row13, new_row14, new_row15, new_row16, new_row17, new_row18, 
              new_row19, new_row20, new_row21)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)



################################ Code Marina

sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","precipitation","mzb.biomass"))
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

################ SERr


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

