setwd("C:/Users/jan/OneDrive/Postdoc HIFMB/Data/Birds")

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

############################################### Housekeeping

data <- read.csv("MARISCO_birds_wadden.csv")

data$StationID <- as.factor(data$StationID)
unique(data$StationID)
levels(data$StationID)

data2 <- filter(data, StationID %in% c("Manslagter Nacken: Manslagter Nacken Sued", 
                                      "Manslagter Nacken: Manslagter Nacken Nord", 
                                      "Leybucht Sued: Leyhoern", 
                                      "Rysumer Nacken: Vorland",
                                      "Pilsum-Manslagt: Pilsum - Manslagter Vorland",
                                      "Leybucht Sued: Hauener Hooge",
                                      "Leybucht Mitte: Mittelplate", 
                                      "Leybucht Mitte: Leybuchtpolder Heller", 
                                      "Leybucht Nord: Buscher Heller", 
                                      "Borkum Sued: Vorland Hopp", 
                                      "Borkum Sued: Binnendeichsflaechen mit Tueskendoersee",
                                      "Borkum Nordwest: Borkum Nordwest",
                                      "Borkum Ostplate: Borkum Ostplate",
                                      "Borkum Sued: Suedstrand - Woldeduenen",
                                      "Norddeich Ost: Vorland",
                                      "Hilgenriedersiel: Sommerpolder",
                                      "Hilgenriedersiel bis Nessmersiel: Sommerpolder",
                                      "Hilgenriedersiel bis Nessmersiel: Vorland",
                                      "Hilgenriedersiel: Vorland",
                                      "Norderney Hafen: Westkopf",
                                      "Norderney Hafen: Hafenbucht",
                                      "Norderney West: Südstrandpolder",
                                      "Norderney Nordwest: Nordbad",
                                      "Norderney Nordwest: Meiereiwiesen",
                                      "Norderney West: Vorland",
                                      "Norderney West: Grohde Polder",
                                      "Norderney Ost: Ostheller - Ostbake",
                                      "Norderney Ost: Nordoststrand",
                                      "Norderney Ost: Nordstrand Mitte",
                                      "Langeoog West: Langeoog West",
                                      "Langeoog West: Ort",
                                      "Langeoog Sued: Langeoog Sued",
                                      "Langeoog Nord: Duenen und Schlopp",
                                      "Langeoog Ost: Vogelkolonie",
                                      "Langeoog Ost: Sommerpolder mit Ostheller",
                                      "Spiekeroog West: Spiekeroog West",
                                      "Spiekeroog West: Spiekeroog Nordweststrand",
                                      "Spiekeroog Ost: Nordoststrand",
                                      "Spiekeroog Ost: Spiekeroog Ost",
                                      "Minsener Oog: Suedwatt",
                                      "Minsener Oog: Nordwatt",
                                      "Minsener Oog: Insel",
                                      "Mellum: Sandbank",
                                      "Mellum: Mellum",
                                      "Robbenplate Eversand: Eversand",
                                      "Dorumer Neufeld Süd: Dorumer Neufeld Süd",
                                      "Wremen Nord: Wremen Nord",
                                      "Langwarder Deich: Vorland",
                                      "Ruhwarden / Tossens: Vorland",
                                      "Burhave: Fedderwardersiel - Burhaversiel",
                                      "Burhave: Burhaversiel - Waddenserdeich",
                                      "Langwarder Deich: Langwarder und Fedderwarder-Groden",
                                      "Jadebusen West: Vorland Nord",
                                      "Jadebusen West: Vorland Sued",
                                      "Jadebusen West: Binnendeichsflaechen",
                                      "Jadebusen Suedost: Vorland Diekmannshausen",
                                      "Jadebusen Suedost: Suederkleihoerne",
                                      "Jadebusen Nord: Jadebusen Nord",
                                      "Jadebusen Ost: Vorland",
                                      "Eckwarderhoerne: Eckwarderhoerne - Iffens",
                                      "Langluetjen: Langluetjen - Blexen",
                                      "Langluetjen: Tettenserhoerne - Langluetjen",
                                      "Waddenserdeich: Waddenserdeich",
                                      "Bremerhaven Sued / Luneplate Ost: Fischereihafen + Außendeich",
                                      "Tegeler Plate / Luneplate West: Tegelerplate / Neues Pfand",
                                      "Bremerhaven Sued / Luneplate Ost: Fischereihafen + Aussendeich",
                                      "Spieka Sued: Vorland",
                                      "Spieka Sued: Sommerpolder",
                                      "Spieka Nord: Vorland",
                                      "Spieka Nord: Sommerpolder",
                                      "Cuxhaven West: Sahlenburg - Arensch Nord"
                                 ))

unique(data2$StationID)
levels(data2$StationID)

rysumer_nacken <- filter(data, StationID %in% c("Manslagter Nacken: Manslagter Nacken Sued",
                                                "Rysumer Nacken: Vorland"))
rysumer_nacken$date<-dmy(rysumer_nacken$date)
rysumer_nacken<-ddply(rysumer_nacken,.(date, species), 
       colwise(mean,  .(abundance, lat, long)), na.rm=T)
rysumer_nacken$lat <- 53.39378
rysumer_nacken$long <- 7.00923

pilsum_manslagt <- filter(data, StationID %in% c("Manslagter Nacken: Manslagter Nacken Nord",
                                                 "Pilsum-Manslagt: Pilsum - Manslagter Vorland"))
pilsum_manslagt$date<-dmy(pilsum_manslagt$date)
pilsum_manslagt<-ddply(pilsum_manslagt,.(date, species), 
                      colwise(mean,  .(abundance, lat, long)), na.rm=T)
pilsum_manslagt$lat <- 53.48254
pilsum_manslagt$long <- 7.02771

leybucht <- filter(data, StationID %in% c("Leybucht Sued: Leyhoern", 
                                       "Leybucht Sued: Hauener Hooge", 
                                       "Leybucht Mitte: Mittelplate", 
                                       "Leybucht Mitte: Leybuchtpolder Heller", 
                                       "Leybucht Nord: Buscher Heller"))
leybucht$date<-dmy(leybucht$date)
leybucht<-ddply(leybucht,.(date, species), 
                       colwise(mean,  .(abundance, lat, long)), na.rm=T)
leybucht$lat <- 53.52555
leybucht$long <- 7.06468

borkum <- filter(data, StationID %in% c("Borkum Sued: Vorland Hopp", 
                                     "Borkum Sued: Binnendeichsflaechen mit Tueskendoersee", 
                                     "Borkum Nordwest: Borkum Nordwest", 
                                     "Borkum Ostplate: Borkum Ostplate", 
                                     "Borkum Sued: Suedstrand - Woldeduenen"))
borkum$date<-dmy(borkum$date)
borkum<-ddply(borkum,.(date, species), 
                colwise(mean,  .(abundance, lat, long)), na.rm=T)
borkum$lat <- 53.56597
borkum$long <- 6.70967

hilgenriedersiel <- filter(data, StationID %in% c("Norddeich Ost: Vorland",
                                               "Hilgenriedersiel: Sommerpolder",
                                               "Hilgenriedersiel bis Nessmersiel: Sommerpolder",
                                               "Hilgenriedersiel bis Nessmersiel: Vorland",
                                               "Hilgenriedersiel: Vorland"))
hilgenriedersiel$date<-dmy(hilgenriedersiel$date)
hilgenriedersiel<-ddply(hilgenriedersiel,.(date, species), 
              colwise(mean,  .(abundance, lat, long)), na.rm=T)
hilgenriedersiel$lat <- 53.66897
hilgenriedersiel$long <- 7.24654

norderney <- filter(data, StationID %in% c("Norderney Hafen: Westkopf",
                                         "Norderney Hafen: Hafenbucht",
                                         "Norderney West: Südstrandpolder",
                                         "Norderney Nordwest: Nordbad",
                                         "Norderney Nordwest: Meiereiwiesen",
                                         "Norderney West: Vorland",
                                         "Norderney West: Grohde Polder",
                                         "Norderney Ost: Ostheller - Ostbake",
                                         "Norderney Ost: Nordoststrand",
                                         "Norderney Ost: Nordstrand Mitte"))
norderney$date<-dmy(norderney$date)
norderney<-ddply(norderney,.(date, species), 
                        colwise(mean,  .(abundance, lat, long)), na.rm=T)
norderney$lat <- 53.72005
norderney$long <- 7.25419

langeoog <- filter(data, StationID %in% c("Langeoog West: Langeoog West",
                                          "Langeoog West: Ort",
                                          "Langeoog Sued: Langeoog Sued",
                                          "Langeoog Nord: Duenen und Schlopp",
                                          "Langeoog Ost: Vogelkolonie",
                                          "Langeoog Ost: Sommerpolder mit Ostheller"))
langeoog$date<-dmy(langeoog$date)
langeoog<-ddply(langeoog,.(date, species), 
                 colwise(mean,  .(abundance, lat, long)), na.rm=T)
langeoog$lat <- 53.74581
langeoog$long <- 7.59503

spiekeroog <- filter(data, StationID %in% c("Spiekeroog West: Spiekeroog West",
                                          "Spiekeroog West: Spiekeroog Nordweststrand",
                                          "Spiekeroog Ost: Nordoststrand",
                                          "Spiekeroog Ost: Spiekeroog Ost"))
spiekeroog$date<-dmy(spiekeroog$date)
spiekeroog<-ddply(spiekeroog,.(date, species), 
                colwise(mean,  .(abundance, lat, long)), na.rm=T)
spiekeroog$lat <- 53.76387
spiekeroog$long <- 7.76623

minsener_oog <- filter(data, StationID %in% c("Minsener Oog: Suedwatt",
                                            "Minsener Oog: Nordwatt",
                                            "Minsener Oog: Insel"))
minsener_oog$date<-dmy(minsener_oog$date)
minsener_oog<-ddply(minsener_oog,.(date, species), 
                  colwise(mean,  .(abundance, lat, long)), na.rm=T)
minsener_oog$lat <- 53.75257
minsener_oog$long <- 8.01525

mellum <- filter(data, StationID %in% c("Mellum: Sandbank",
                                      "Mellum: Mellum"))
mellum$date<-dmy(mellum$date)
mellum<-ddply(mellum,.(date, species), 
                    colwise(mean,  .(abundance, lat, long)), na.rm=T)
mellum$lat  <- 53.7216
mellum$long <- 8.15202

eversand <- filter(data, StationID %in% c("Robbenplate Eversand: Eversand",
                                          "Dorumer Neufeld Süd: Dorumer Neufeld Süd",
                                          "Wremen Nord: Wremen Nord"))
eversand$date<-dmy(eversand$date)
eversand<-ddply(eversand,.(date, species), 
              colwise(mean,  .(abundance, lat, long)), na.rm=T)
eversand$lat <- 53.70441
eversand$long <- 8.36752

langenwarder_deich <- filter(data, StationID %in% c( "Langwarder Deich: Vorland",
                                           "Ruhwarden / Tossens: Vorland",
                                           "Burhave: Fedderwardersiel - Burhaversiel",
                                           "Burhave: Burhaversiel - Waddenserdeich",
                                           "Langwarder Deich: Langwarder und Fedderwarder-Groden"))
langenwarder_deich$date<-dmy(langenwarder_deich$date)
langenwarder_deich<-ddply(langenwarder_deich,.(date, species), 
                colwise(mean,  .(abundance, lat, long)), na.rm=T)
langenwarder_deich$lat <- 53.61271
langenwarder_deich$long <- 8.31263

jadebusen_west <- filter(data, StationID %in% c("Jadebusen West: Vorland Nord",
                                           "Jadebusen West: Vorland Sued",
                                           "Jadebusen West: Binnendeichsflaechen"))
jadebusen_west$date<-dmy(jadebusen_west$date)
jadebusen_west<-ddply(jadebusen_west,.(date, species), 
                          colwise(mean,  .(abundance, lat, long)), na.rm=T)
jadebusen_west$lat <- 53.49431
jadebusen_west$long <- 8.0611

jadebusen_ost <- filter(data, StationID %in% c("Jadebusen Suedost: Vorland Diekmannshausen",
                                                "Jadebusen Suedost: Suederkleihoerne",
                                                "Jadebusen Nord: Jadebusen Nord",
                                                "Jadebusen Ost: Vorland",
                                                "Eckwarderhoerne: Eckwarderhoerne - Iffens"))
jadebusen_ost$date<-dmy(jadebusen_ost$date)
jadebusen_ost<-ddply(jadebusen_ost,.(date, species), 
                      colwise(mean,  .(abundance, lat, long)), na.rm=T)
jadebusen_ost$lat <- 53.47611
jadebusen_ost$long <- 8.31328

langluetjen <- filter(data, StationID %in% c("Langluetjen: Langluetjen - Blexen",
                                              "Langluetjen: Tettenserhoerne - Langluetjen",
                                              "Waddenserdeich: Waddenserdeich"))
langluetjen$date<-dmy(langluetjen$date)
langluetjen<-ddply(langluetjen,.(date, species), 
                     colwise(mean,  .(abundance, lat, long)), na.rm=T)
langluetjen$lat <- 53.54713
langluetjen$long <- 8.52331

bremerhaven <- filter(data, StationID %in% c("Bremerhaven Sued / Luneplate Ost: Fischereihafen + Außendeich",
                                             "Tegeler Plate / Luneplate West: Tegelerplate / Neues Pfand",
                                             "Bremerhaven Sued / Luneplate Ost: Fischereihafen + Aussendeich"))
bremerhaven$date<-dmy(bremerhaven$date)
bremerhaven<-ddply(bremerhaven,.(date, species), 
                   colwise(mean,  .(abundance, lat, long)), na.rm=T)
bremerhaven$lat <- 53.50341
bremerhaven$long <- 8.5402

spieka <- filter(data, StationID %in% c("Spieka Sued: Vorland",
                                        "Spieka Sued: Sommerpolder",
                                        "Spieka Nord: Vorland",
                                        "Spieka Nord: Sommerpolder",
                                        "Cuxhaven West: Sahlenburg - Arensch Nord"))
spieka$date<-dmy(spieka$date)
spieka<-ddply(spieka,.(date, species), 
                   colwise(mean,  .(abundance, lat, long)), na.rm=T)
spieka$lat <- 53.8127
spieka$long <- 8.56423

data2 <- bind_rows(list(rysumer_nacken=rysumer_nacken, pilsum_manslagt=pilsum_manslagt, 
                        leybucht=leybucht, borkum=borkum, hilgenriedersiel=hilgenriedersiel,
                        norderney=norderney, langeoog=langeoog, spiekeroog=spiekeroog,
                        minsener_oog=minsener_oog, mellum=mellum, eversand=eversand,
                        langenwarder_deich=langenwarder_deich, jadebusen_west=jadebusen_west,
                        jadebusen_ost=jadebusen_ost, langluetjen=langluetjen, 
                        bremerhaven=bremerhaven, spieka=spieka), .id = "StationID")
data2$year <- year(data2$date)


###################################################### alpha-diversity


summary(data2)
data2<-ddply(data2,.(StationID, year, species), 
             colwise(sum,  .(abundance)), na.rm=T)
unique(data2$species)
unique(data2$year)
unique(data2$StationID)

data_wideBM<-dcast(data2, StationID + year ~ species, value.var="abundance")

data_wideBM[,3:266][is.na(data_wideBM[,3:266])]<-0

data <- data_wideBM
length(unique(data$StationID))
max(data$year)-min(data$year)+1
nrow(data)
ncol(data)-2
data$S<-specnumber(data[,3:ncol(data)])
data$ENS<-diversity(data[,3:ncol(data)], "invsimpson")

data<- ddply(data, .(StationID, year), 
                      summarise, 
                      S = mean(S), ENS = mean(ENS))

#write.csv(data, file = "MARISCO_birds_wadden_alpha-div.csv")

###################################################### gamma-diversity

summary(data2)
data2<-ddply(data2,.(StationID, year, species), 
            colwise(mean,  .(abundance)), na.rm=T)
unique(data2$species)
unique(data2$year)
unique(data2$StationID)

data_wideBM<-dcast(data2, StationID + year ~ species, value.var="abundance")

data_wideBM[,3:266][is.na(data_wideBM[,3:266])]<-0

data <- data_wideBM
length(unique(data$StationID))
max(data$year)-min(data$year)+1
nrow(data)
ncol(data)-2
data$S<-specnumber(data[,3:ncol(data)])
data$D<-diversity(data[,3:ncol(data)], "invsimpson")

#write.csv(data, file = "MARISCO_birds_wadden_gamma-div.csv")

########################################## Env-data: MZB biomass

data <- read.csv("MARISCO_mzb_wadden_raw.csv")

rysumer_nacken <- filter(data, StationID %in% c("EmDo_MZB_2",
                                                "EmDo_MZB_5"))
rysumer_nacken$date<-dmy(rysumer_nacken$date)
rysumer_nacken<-ddply(rysumer_nacken,.(date, StationID), 
                      colwise(sum,  .(biomass, lat, long)), na.rm=T)
rysumer_nacken$lat <- 53.39378
rysumer_nacken$long <- 7.00923

pilsum_manslagt <- filter(data, StationID %in% c("Bork_MZB_4",
                                                 "Leyb_MZB_1"))
pilsum_manslagt$date<-dmy(pilsum_manslagt$date)
pilsum_manslagt<-ddply(pilsum_manslagt,.(date, StationID), 
                       colwise(sum,  .(biomass, lat, long)), na.rm=T)
pilsum_manslagt$lat <- 53.48254
pilsum_manslagt$long <- 7.02771

leybucht <- filter(data, StationID %in% c("Leyb_MZB_2", 
                                          "Leyb_MZB_3", 
                                          "Leyb_MZB_7", 
                                          "Leyb_MZB_6", 
                                          "Leyb_MZB_4",
                                          "Leyb_MZB_5"))
leybucht$date<-dmy(leybucht$date)
leybucht<-ddply(leybucht,.(date, StationID), 
                colwise(sum,  .(biomass, lat, long)), na.rm=T)
leybucht$lat <- 53.52555
leybucht$long <- 7.06468

borkum <- filter(data, StationID %in% c("Bork_MZB_8"))

borkum$date<-dmy(borkum$date)
borkum<-ddply(borkum,.(date, StationID), 
              colwise(sum,  .(biomass, lat, long)), na.rm=T)
borkum$lat <- 53.56597
borkum$long <- 6.70967

hilgenriedersiel <- filter(data, StationID %in% c("Nney_MZB_8",
                                                  "Nney_MZB_7",
                                                  "Nney_MZB_6"))
hilgenriedersiel$date<-dmy(hilgenriedersiel$date)
hilgenriedersiel<-ddply(hilgenriedersiel,.(date, StationID), 
                        colwise(sum,  .(biomass, lat, long)), na.rm=T)
hilgenriedersiel$lat <- 53.66897
hilgenriedersiel$long <- 7.24654

norderney <- filter(data, StationID %in% c("Nney_MZB_5",
                                           "Nney_MZB_4",
                                           "Nney_MZB_3",
                                           "Nney_MZB_2",
                                           "Nney_MZB_1",
                                           "Nney_MZB_9"))
norderney$date<-dmy(norderney$date)
norderney<-ddply(norderney,.(date, StationID), 
                 colwise(sum,  .(biomass, lat, long)), na.rm=T)
norderney$lat <- 53.72005
norderney$long <- 7.25419

langeoog <- filter(data, StationID %in% c("Spog_MZB_5"))

langeoog$date<-dmy(langeoog$date)
langeoog<-ddply(langeoog,.(date, StationID), 
                colwise(sum,  .(biomass, lat, long)), na.rm=T)
langeoog$lat <- 53.74581
langeoog$long <- 7.59503

spiekeroog <- filter(data, StationID %in% c("Spog_MZB_2"))

spiekeroog$date<-dmy(spiekeroog$date)
spiekeroog<-ddply(spiekeroog,.(date, StationID), 
                  colwise(sum,  .(biomass, lat, long)), na.rm=T)
spiekeroog$lat <- 53.76387
spiekeroog$long <- 7.76623

minsener_oog <- filter(data, StationID %in% c("Jade_MZB_1"))

minsener_oog$date<-dmy(minsener_oog$date)
minsener_oog<-ddply(minsener_oog,.(date, StationID), 
                    colwise(sum,  .(biomass, lat, long)), na.rm=T)
minsener_oog$lat <- 53.75257
minsener_oog$long <- 8.01525

mellum <- filter(data, StationID %in% c("JaBu_MZB_12"))

mellum$date<-dmy(mellum$date)
mellum<-ddply(mellum,.(date, StationID), 
              colwise(sum,  .(biomass, lat, long)), na.rm=T)
mellum$lat  <- 53.7216
mellum$long <- 8.15202

eversand <- filter(data, StationID %in% c("AuWe_MZB_1"))

eversand$date<-dmy(eversand$date)
eversand<-ddply(eversand,.(date, StationID), 
                colwise(sum,  .(biomass, lat, long)), na.rm=T)
eversand$lat <- 53.70441
eversand$long <- 8.36752

langenwarder_deich <- filter(data, StationID %in% c( "AuWe_MZB_3"))

langenwarder_deich$date<-dmy(langenwarder_deich$date)
langenwarder_deich<-ddply(langenwarder_deich,.(date, StationID), 
                          colwise(sum,  .(biomass, lat, long)), na.rm=T)
langenwarder_deich$lat <- 53.61271
langenwarder_deich$long <- 8.31263

jadebusen_west <- filter(data, StationID %in% c("JaBu_MZB_9"))

jadebusen_west$date<-dmy(jadebusen_west$date)
jadebusen_west<-ddply(jadebusen_west,.(date, StationID), 
                      colwise(sum,  .(biomass, lat, long)), na.rm=T)
jadebusen_west$lat <- 53.49431
jadebusen_west$long <- 8.0611

jadebusen_ost <- filter(data, StationID %in% c("JaBu_MZB_8"))

jadebusen_ost$date<-dmy(jadebusen_ost$date)
jadebusen_ost<-ddply(jadebusen_ost,.(date, StationID), 
                     colwise(sum,  .(biomass, lat, long)), na.rm=T)
jadebusen_ost$lat <- 53.47611
jadebusen_ost$long <- 8.31328

langluetjen <- filter(data, StationID %in% c("WeMu_MZB_1",
                                             "WeMu_MZB_2",
                                             "WeMu_MZB_3",
                                             "WeMu_MZB_4"))

langluetjen$date<-dmy(langluetjen$date)
langluetjen<-ddply(langluetjen,.(date, StationID), 
                   colwise(sum,  .(biomass, lat, long)), na.rm=T)
langluetjen$lat <- 53.54713
langluetjen$long <- 8.52331

bremerhaven <- filter(data, StationID %in% c("WeMu_MZB_7",
                                             "WeMu_MZB_6",
                                             "WeMu_MZB_5"))
bremerhaven$date<-dmy(bremerhaven$date)
bremerhaven<-ddply(bremerhaven,.(date, StationID), 
                   colwise(sum,  .(biomass, lat, long)), na.rm=T)
bremerhaven$lat <- 53.50341
bremerhaven$long <- 8.5402

spieka <- filter(data, StationID %in% c("WuKu_MZB_6"))

spieka$date<-dmy(spieka$date)
spieka<-ddply(spieka,.(date, StationID), 
              colwise(sum,  .(biomass, lat, long)), na.rm=T)
spieka$lat <- 53.8127
spieka$long <- 8.56423

data2 <- bind_rows(list(rysumer_nacken=rysumer_nacken, pilsum_manslagt=pilsum_manslagt, 
                        leybucht=leybucht, borkum=borkum, hilgenriedersiel=hilgenriedersiel,
                        norderney=norderney, langeoog=langeoog, spiekeroog=spiekeroog,
                        minsener_oog=minsener_oog, mellum=mellum, eversand=eversand,
                        langenwarder_deich=langenwarder_deich, jadebusen_west=jadebusen_west,
                        jadebusen_ost=jadebusen_ost, langluetjen=langluetjen, 
summary(data2)
data2<-ddply(data2,.(StationID, date), 
             colwise(sum,  .(biomass)), na.rm=T)
data2$year <- year(data2$date)
data2<-ddply(data2,.(StationID, year), 
             colwise(mean,  .(biomass)), na.rm=T)

unique(data2$year)
unique(data2$StationID)

#write.csv(data2, file = "MARISCO_birds_wadden_MZBbiomass.csv")

data <- read.csv("MARISCO_birds_wadden_SEM.csv")
data2 <- read.csv("MARISCO_birds_wadden_MZBbiomass.csv")

bird_combo <- merge(data, data2, by=c("StationID", "year"))

#write.csv(bird_combo, file = "MARISCO_birds_wadden_SEM.csv")

####################################################### SEM Analyis

setwd("C:/Users/jan/OneDrive/Postdoc HIFMB/Data/Birds")

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


####################################################### GDM Analyis

data2 <- bind_rows(list(rysumer_nacken=rysumer_nacken, pilsum_manslagt=pilsum_manslagt, 
                        leybucht=leybucht, borkum=borkum, hilgenriedersiel=hilgenriedersiel,
                        norderney=norderney, langeoog=langeoog, spiekeroog=spiekeroog,
                        minsener_oog=minsener_oog, mellum=mellum, eversand=eversand,
                        langenwarder_deich=langenwarder_deich, jadebusen_west=jadebusen_west,
                        jadebusen_ost=jadebusen_ost, langluetjen=langluetjen, 
                        bremerhaven=bremerhaven, spieka=spieka), .id = "StationID")
data2$year <- year(data2$date)

data2<-ddply(data2,.(StationID, year, species), 
             colwise(mean,  .(abundance)), na.rm=T)
unique(data2$species)
unique(data2$year)
unique(data2$StationID)

################################################## Rysumer Nacken

data <- filter(data2, StationID %in% c("rysumer_nacken"))

data_wide<-dcast(data, year ~ species, value.var="abundance")

data_wide[,2:82][is.na(data_wide[,2:82])]<-0

#write.csv(data_wide, file = "rysumer_nacken_wide.csv")

################################################## Pilsum Manslagt

data <- filter(data2, StationID %in% c("pilsum_manslagt"))

data_wide<-dcast(data, year ~ species, value.var="abundance")

data_wide[,2:79][is.na(data_wide[,2:79])]<-0

#write.csv(data_wide, file = "pilsum_manslagt_wide.csv")

################################################## Leybucht

data <- filter(data2, StationID %in% c("leybucht"))

data_wide<-dcast(data, year ~ species, value.var="abundance")

data_wide[,2:151][is.na(data_wide[,2:151])]<-0

#write.csv(data_wide, file = "leybucht_wide.csv")

################################################## Borkum

data <- filter(data2, StationID %in% c("borkum"))

data_wide<-dcast(data, year ~ species, value.var="abundance")

data_wide[,2:157][is.na(data_wide[,2:157])]<-0

#write.csv(data_wide, file = "borkum_wide.csv")

################################################## Hilgenriedersiel

data <- filter(data2, StationID %in% c("hilgenriedersiel"))

data_wide<-dcast(data, year ~ species, value.var="abundance")

data_wide[,2:154][is.na(data_wide[,2:154])]<-0

#write.csv(data_wide, file = "hilgenriedersiel_wide.csv")

################################################## Norderney

data <- filter(data2, StationID %in% c("norderney"))

data_wide<-dcast(data, year ~ species, value.var="abundance")

data_wide[,2:162][is.na(data_wide[,2:162])]<-0

#write.csv(data_wide, file = "norderney_wide.csv")

################################################## Langeoog

data <- filter(data2, StationID %in% c("langeoog"))

data_wide<-dcast(data, year ~ species, value.var="abundance")

data_wide[,2:161][is.na(data_wide[,2:161])]<-0

#write.csv(data_wide, file = "langeoog_wide.csv")

################################################## Spiekeroog

data <- filter(data2, StationID %in% c("spiekeroog"))

data_wide<-dcast(data, year ~ species, value.var="abundance")

data_wide[,2:134][is.na(data_wide[,2:134])]<-0

#write.csv(data_wide, file = "spiekeroog_wide.csv")

################################################## Minsener Oog

data <- filter(data2, StationID %in% c("minsener_oog"))

data_wide<-dcast(data, year ~ species, value.var="abundance")

data_wide[,2:126][is.na(data_wide[,2:126])]<-0

#write.csv(data_wide, file = "minsener_oog_wide.csv")

################################################## Mellum

data <- filter(data2, StationID %in% c("mellum"))

data_wide<-dcast(data, year ~ species, value.var="abundance")

data_wide[,2:140][is.na(data_wide[,2:140])]<-0

#write.csv(data_wide, file = "mellum_wide.csv")

################################################## Jadebusen West

data <- filter(data2, StationID %in% c("jadebusen_west"))

data_wide<-dcast(data, year ~ species, value.var="abundance")

data_wide[,2:186][is.na(data_wide[,2:186])]<-0

#write.csv(data_wide, file = "jadebusen_west_wide.csv")

################################################## Jadebusen Ost

data <- filter(data2, StationID %in% c("jadebusen_ost"))

data_wide<-dcast(data, year ~ species, value.var="abundance")

data_wide[,2:109][is.na(data_wide[,2:109])]<-0

#write.csv(data_wide, file = "jadebusen_ost_wide.csv")

################################################## Jadebusen Ost

data <- filter(data2, StationID %in% c("langenwarder_deich"))

data_wide<-dcast(data, year ~ species, value.var="abundance")

data_wide[,2:93][is.na(data_wide[,2:93])]<-0

#write.csv(data_wide, file = "langenwarder_deich_wide.csv")

################################################## Langlütjen

data <- filter(data2, StationID %in% c("langluetjen"))

data_wide<-dcast(data, year ~ species, value.var="abundance")

data_wide[,2:99][is.na(data_wide[,2:99])]<-0

#write.csv(data_wide, file = "langluetjen_wide.csv")

################################################## Bremerhaven

data <- filter(data2, StationID %in% c("bremerhaven"))

data_wide<-dcast(data, year ~ species, value.var="abundance")

data_wide[,2:96][is.na(data_wide[,2:96])]<-0

#write.csv(data_wide, file = "bremerhaven_wide.csv")

################################################## Spieka

data <- filter(data2, StationID %in% c("spieka"))

data_wide<-dcast(data, year ~ species, value.var="abundance")

data_wide[,2:95][is.na(data_wide[,2:95])]<-0

#write.csv(data_wide, file = "spieka_wide.csv")

################################################## Eversand

data <- filter(data2, StationID %in% c("eversand"))

data_wide<-dcast(data, year ~ species, value.var="abundance")

data_wide[,2:60][is.na(data_wide[,2:60])]<-0

#write.csv(data_wide, file = "eversand_wide.csv")

################################################# Selected sites

data <- filter(data2, StationID %in% c("rysumer_nacken", "borkum", "leybucht", "hilgenriedersiel", "norderney", 
                                       "pilsum_manslagt", "bremerhaven", "jadebusen_ost", "jadebusen_west",
                                       "langenwarder_deich", "langeoog", "langluetjen"))

data_wide<-dcast(data, year + StationID ~ species, value.var="abundance")

data_wide[,3:261][is.na(data_wide[,3:261])]<-0

#write.csv(data_wide, file = "selected_sites_wide.csv")

data <- read.csv("selected_sites_wide.csv")

data <- summarise_each(group_by(data, year), funs(mean(., na.rm = TRUE)))

write.csv(data, file = "select_wide.csv")


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

