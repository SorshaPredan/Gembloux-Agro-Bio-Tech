# PC Gembloux
setwd("C:/Users/sorsha/Desktop/Analisi Dendro 2026")
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

# SPEI SUMMER
spei_data <- data.frame(
  year = PrecSites$year,
  Jan = PrecSites$Jan - TempSites$Jan,
  Feb = PrecSites$Feb - TempSites$Feb,
  Mar = PrecSites$Mar - TempSites$Mar,
  Apr = PrecSites$Apr - TempSites$Apr,
  May = PrecSites$May - TempSites$May,
  Jun = PrecSites$Jun - TempSites$Jun,
  Jul = PrecSites$Jul - TempSites$Jul,
  Aug = PrecSites$Aug - TempSites$Aug,
  Sep = PrecSites$Sep - TempSites$Sep,
  Oct = PrecSites$Oct - TempSites$Oct,
  Nov = PrecSites$Nov - TempSites$Nov,
  Dec = PrecSites$Dec - TempSites$Dec
)
head(spei_data)
spei_summer_data <- data.frame(
  year = spei_data$year,
  summer = rowMeans(
    spei_data[,c("Jun","Jul","Aug")],
    na.rm = TRUE
  )
)
head(spei_summer_data)
spei_summer <- spei(
  ts(
    spei_summer_data$summer,
    start = min(spei_summer_data$year),
    frequency = 1
  ),
  scale = 1
)
# Lo SPEI normalmente:
# 0 = normale
# negativo = siccità
# -1 = siccità moderata
# -1.5 = siccità severa
# -2 = siccità estrema
drought_years <- data.frame(
  year = spei_summer_data$year,
  spei = as.numeric(spei_summer$fitted)
)
head(drought_years)
extreme_drought <- drought_years[
  drought_years$spei <= -1.5,
]
extreme_drought
nrow(extreme_drought)

# SEA (Superposed Epoch Analysis)
extreme_drought$year %in% as.numeric(rownames(EWDChron))
drought_events <- c(
  1922,
  1927,
  1953,
  1965,
  1972,
  1980,
  1987,
  2007,
  2012
)
drought_events

class(EWD) <- c("rwl","data.frame")
class(LWD) <- c("rwl","data.frame")
class(MXD) <- c("rwl","data.frame")
EWD_SEA <- sea(
  EWD,
  key = drought_events,
  lag = 3,
  resample = 1000
)
plot(EWD_SEA)

LWD_SEA <- sea(
  LWD,
  key = drought_events,
  lag = 3,
  resample = 1000
)
plot(LWD_SEA)

MXD_SEA <- sea(
  MXD,
  key = drought_events,
  lag = 3,
  resample = 1000
)
plot(MXD_SEA)
range(as.numeric(rownames(EWDChron)))
range(PrecSites$year)
range(TempSites$year)
