
# PC Sorsha
setwd("C:/Users/user/Desktop/TIROCINIO FINALE/ANALISI DENDRO/pourSorsha")
remove(list = ls())

library("dplR")
library("tidyverse")
library("dplyr")
library("treeclim")
library("ggplot2")
library("SPEI")
library("pointRes")
library("bootRes") # Watch out, the DCC command is the same as in treeclim, but different arguments!
library("corrplot")
library("Rcpp")
library(data.table)
library(stats)
library(knitr)
library(graphics)
library(utils)
#### Spatial Correlation ####

library(ncdf4)
library(fields)
library(Hmisc)
library(mapdata)
library(readxl)
# Read in MXD
install.packages("remotes")
remotes::install_github("AllanBuras/dendRolAB")

# DATI CLIMATICI
PrecSites<-read.table("Prec Site.txt", header = TRUE)
#temp <- read.table("tempBISHOP.txt", header = FALSE)
TempSites<-read.table("Temp Site.txt", header=TRUE)

# CREAZIONE FILE TUCSON CORRETTO PER COFECHA
## IMPORT DEL FILE ORIGINALE
TRW <- read_excel(
  file.choose(),
  col_names = FALSE
)
TRW <- as.data.frame(TRW)
head(TRW)
dim(TRW)
# CONTROLLO STRUTTURA
# La prima riga contiene i nomi dei campioni
# La prima colonna contiene gli anni (YEARS)
campioni <- as.character(TRW[1, -1])
anni <- TRW[-1, 1]

# CREAZIONE MATRICE RWL
TRW <- TRW[-1, -1]
names(TRW) <- campioni
TRW[] <- lapply(TRW, as.numeric)
rownames(TRW) <- anni
TRW[TRW == 999] <- NA
class(TRW) <- c("rwl", "data.frame")
head(TRW)
dim(TRW)

# ESPORTAZIONE FILE TUCSON PER COFECHA
write.rwl(
  TRW,
  "CAMPIONI.rwl",
  format = "tucson",
  long.names = TRUE
)

# ANALISI DENDROCRONOLOGICA
# ricarico il file Tucson pulito
# PC Gembloux
# NO PC Sorsha                     
TRW <- read.rwl("campioni.rwl")
head(TRW)
# Riprendere da qui                     
# Controllo qualità
rwi.stats(TRW)
rwi.stats.running(TRW)
corr.rwl.seg(TRW)

library(openxlsx)
write.xlsx(TRW, "C:/Users/user/Desktop/TIROCINIO FINALE/ANALISI DENDRO/pourSorsha/TRW.xlsx", rowNames = TRUE)
ls()

# BAI (Basal Area Increment = incrementO dell’area basale)
## gli alberi malati crescono meno?
BeechBAI <- bai.in(TRW)
# Grafico serie grezze
matplot(
  as.numeric(rownames(TRW)),
  TRW,
  type="l",
  xlab="Year",
  ylab="Tree ring width"
)                      
# Standardizzazione
TRWdetrend<-detrend(TRW, method = "Spline", nyrs = 30)
# Cronologia                      
BeechChron<-chron(TRWdetrend,prefix = "AVG", biweight = TRUE, prewhiten = FALSE)
plot.crn(BeechChron)                      
range(time(BeechChron))

# SITES
# Precipitation
class(PrecSites)
head(PrecSites)
dim(PrecSites)
colnames(PrecSites) <- c("year",
                         "Jan","Feb","Mar","Apr","May","Jun",
                         "Jul","Aug","Sep","Oct","Nov","Dec")
PrecSites[,2:13] <- lapply(
  PrecSites[,2:13],
  function(x) as.numeric(gsub("\\.", "", x))
)
PrecSites[,2:13] <- PrecSites[,2:13] / 1000000
summary(PrecSites[,2:13])
duplicati <- PrecSites$year[duplicated(PrecSites$year)]
PrecSites[PrecSites$year %in% duplicati, ]
PrecSites <- PrecSites[!duplicated(PrecSites$year), ]
plot(dcc(BeechChron, PrecSites, selection = -6:9, method = "correlation",
           timespan = c(1930,1990), var_names = "precipitation", boot = "std"))

# Temperature
class(TempSites)
head(TempSites)
dim(TempSites)
colnames(TempSites) <- c("year",
                         "Jan","Feb","Mar","Apr","May","Jun",
                         "Jul","Aug","Sep","Oct","Nov","Dec")
TempSites[,2:13] <- lapply(
  TempSites[,2:13],
  function(x) as.numeric(gsub("\\.", "", x))
)
TempSites[,2:13] <- TempSites[,2:13] / 1000000   
summary(TempSites[,2:13])
duplicati_temp <- TempSites$year[duplicated(TempSites$year)]
TempSites[TempSites$year %in% duplicati_temp, ]
TempSites <- TempSites[!duplicated(TempSites$year), ]
plot(dcc(BeechChron, TempSites, selection = -6:9, method = "correlation",
           timespan = c(1930,1990), var_names = "temperature", boot = "std"))

# CORRELAZIONE CLIMA
## CORRELAZIONE TRW - PRECIPITAZIONI
PrecCorr <- dcc(BeechChron, PrecSites, selection = -6:9, method = "correlation",
                  timespan = c(1930,1990), var_names = "precipitation",boot = "std")
## CORRELAZIONE TRW - TEMPERATURA
TempCorr <- dcc(BeechChron,TempSites,selection = -6:9,method = "correlation",
                  timespan = c(1930,1990), var_names = "temperature", boot = "std")
