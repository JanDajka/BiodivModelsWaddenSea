setwd(".../Data/Fish")

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

###################################################################### Standing Diversity Model

data <- read.csv("MARISCO_fish_wadden_SEM.csv")

data$StationID <- as.factor(data$StationID)

#Are the missing values?

levels(data$StationID)
colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

summary(data)

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

boxplot(data$biomass, 
        main = "Biomass")
hist(data$biomass)
dotchart(data$biomass)

par(mfrow = c(1, 2))
plot(x = data$year, 
     y = data$S)
plot(x = data$year, 
     y = data$ENS)
par(mfrow = c(1, 1))
plot(x = data$year, 
     y = data$biomass)

# Collinearity X

names(data)
pairs.panels(data[, c("temperature", "oxygen", "pH", 
                          "S", "ENS", "biomass")])

# Zero inflation
sum(data == 0)

# Scale variables (method for GLMs)
data$temperature <- scale(data$temperature)
data$oxygen <- scale(data$oxygen)
data$pH <- scale(data$pH)

# Recode vars to roughly same scale
#data$biomass <- data$biomass/10000
data$biomass <- log(data$biomass)
data$S <- data$S/10


################################## gls fitting


gls.abu.S <- gls(S ~ temperature + oxygen, data=data,
                 correlation=corGaus(form= ~ lat + long | year))
summary(gls.abu.S)
plot(gls.abu.S)
hist(resid(gls.abu.S,type="pearson"))

gls.abu.ENS <- gls(ENS ~ temperature + pH, data=data,
                   correlation=corGaus(form= ~ long + lat | year))
summary(gls.abu.ENS)
plot(gls.abu.ENS)
hist(resid(gls.abu.ENS,type="pearson"))

gls.abu.biom <- gls(biomass ~ temperature + oxygen + ENS + S, data=data,
                    correlation=corGaus(form= ~ long + lat | year))
summary(gls.abu.biom)
plot(gls.abu.biom)
hist(resid(gls.abu.biom,type="pearson"))


sem.SD.abu.gls1 <- psem(
  gls.abu.S,
  gls.abu.ENS,
  gls.abu.biom,
  S %~~% ENS,
  data)
summary(sem.SD.abu.gls1)
plot(sem.SD.abu.gls1)


################################################## GDM Model

##################### SERr

source("turnover.R")

data <- read.csv("Balt_W_1_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 8:N;


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
SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10, new_row11, new_row12, new_row13, new_row14)
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
SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10, new_row11, new_row12, new_row13, new_row14)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)



sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","pH","oxygen"))
pred <- names(data[,env])
d <- data

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

##################### SERa

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

##################### SERr

source("turnover.R")

data <- read.csv("Jade_W_1_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 8:N;


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
SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10, new_row11, new_row12, new_row13, new_row14)
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
SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10, new_row11, new_row12, new_row13, new_row14)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)



sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","pH","oxygen"))
pred <- names(data[,env])
d <- data

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


##################### SERa

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

source("turnover.R")

data <- read.csv("Spog_W_1_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 8:N;


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
SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10, new_row11, new_row12, new_row13, new_row14)
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
SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10, new_row11, new_row12, new_row13, new_row14)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","pH","oxygen"))
pred <- names(data[,env])
d <- data

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

##################### SERa

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

source("turnover.R")

data <- read.csv("Spog_W_2_wide.csv")

colSums(is.na(data))

data <- na.exclude(data) 

dim(data)

#number of rows and columns
M = nrow(data)
N = ncol(data)

#Define columns with species abundance data
SpecColumns = 8:N;


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
SERr <- rbind(SERr, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10, new_row11, new_row12, new_row13, new_row14)
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
SERa <- rbind(SERa, new_row1, new_row2, new_row3, new_row4, new_row5, new_row6, new_row7, new_row8, new_row9, new_row10, new_row11, new_row12, new_row13, new_row14)
SERa2 <- acast(SERa, To ~ From , value.var="SER", na.rm= TRUE)


sp <- SpecColumns
env <- which(colnames(data) %in% c("temperature","pH","oxygen"))
pred <- names(data[,env])
d <- data

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

##################### SERa

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
