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
prec<-read.table("Prec Chimay.txt", header = TRUE)
#temp <- read.table("tempBISHOP.txt", header = FALSE)
temp<-read.table("Temp Uccle.txt", header=TRUE)

PrecSites<-read.table("Prec Site.txt", header = TRUE)
#temp <- read.table("tempBISHOP.txt", header = FALSE)
TempSites<-read.table("Temp Site.txt", header=TRUE)

# CREAZIONE FILE TUCSON CORRETTO PER COFECHA
## IMPORT DEL FILE ORIGINALE: (il tuo campioni.rwl è letto come testo perché non è fixed width)
raw <- read.table(
  "campioni.rwl",
  header = FALSE,
  fill = TRUE
)
# Ricostruzione serie annuali
campioni <- unique(raw$V1)
serie <- list()
for(s in campioni){
  x <- raw[raw$V1 == s, ]
  anni <- c()
  valori <- c()
  for(i in 1:nrow(x)){
    anno_inizio <- x$V2[i]
    valori_riga <- as.numeric(x[i,3:ncol(x)])
    anni <- c(
      anni,
      anno_inizio + 0:(length(valori_riga)-1)
    )
    valori <- c(
      valori,
      valori_riga)
  }
serie[[s]] <- data.frame(
    anno = anni,
    valore = valori)
}

# CONTROLLO ANNI ANOMALI
# individua eventuali ID concatenati con l'anno
anni_min <- sapply(serie, function(x) min(x$anno))
problemi <- names(anni_min[anni_min < 1000])
if(length(problemi) > 0){
  print("Attenzione: possibili campioni con anno errato:")
  print(problemi)
}                
# correzione del caso trovato nel file
if("FSMA134A1860" %in% names(serie)){
  serie[["FSMA134A"]] <- serie[["FSMA134A1860"]]
  serie[["FSMA134A"]]$anno <-
    serie[["FSMA134A"]]$anno + 1786
  serie[["FSMA134A1860"]] <- NULL
}
campioni <- names(serie)

# CREAZIONE OGGETTO RWL
anni <- unlist(lapply(serie, function(x)x$anno))
anni_totali <- min(anni):max(anni)
TRW <- matrix(
  NA,
  nrow = length(anni_totali),
  ncol = length(campioni),
  dimnames=list(
    as.character(anni_totali),
    campioni)
)
for(s in campioni){
  TRW[
    as.character(serie[[s]]$anno),
    s
  ] <-
    serie[[s]]$valore
}

# codici mancanti Tucson
TRW[TRW == 999] <- NA
TRW <- as.data.frame(TRW)
class(TRW) <- c("rwl","data.frame")

# ESPORTAZIONE TUCSON PER COFECHA
write.rwl(
  TRW,
  "CAMPIONI.rwl",
  format = "tucson",
  long.names = TRUE
)

# ANALISI DENDROCRONOLOGICA
# ricarico il file Tucson pulito
TRW <- read.rwl(
  "campioni.rwl"
)
head(TRW)
# Controllo qualità
rwi.stats(TRW)
rwi.stats.running(TRW)
corr.rwl.seg(TRW)

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
                      



# BAI Becch Malade/No Malade
## gli alberi malati presentano una crescita inferiore rispetto agli alberi sani?
HealthyTRW <- read_xlsx("TRWBeechNoMalade.xlsx") 
HealthyTRW <- as.data.frame(HealthyTRW)
DiseasedTRW <- read_xlsx("TRWBeechMalade.xlsx")
DiseasedTRW <- as.data.frame(DiseasedTRW)
HealthyTRW[1] <- NULL
DiseasedTRW[1] <- NULL
row.names(HealthyTRW)<-1819:2025
row.names(DiseasedTRW)<-1819:2025
HealthyBAI <- bai.in(HealthyTRW)
years_HealthyBAI <- as.numeric(row.names(HealthyBAI))
DiseasedBAI <- bai.in(DiseasedTRW)
years_DiseasedBAI <- as.numeric(row.names(DiseasedBAI))
## Media annuale degli alberi sani:
mean_HealthyBAI <- rowMeans(HealthyBAI, na.rm = TRUE)
## Media annuale degli alberi malati:
mean_DiseasedBAI <- rowMeans(DiseasedBAI, na.rm = TRUE)
plot( years_HealthyBAI, mean_HealthyBAI, type="l", col="blue", lwd=2, xlab="Year", ylab="BAI (cm²/anno)" ) 
      lines( years_DiseasedBAI, mean_DiseasedBAI, col="red", lwd=2 ) 
     legend( "topright", legend=c("Beech No Malade", "Beech Malade"), col=c("blue","red"), lwd=2 )
### Differenza di crescita tra i gruppi
difference_BAI <- mean_HealthyBAI - mean_DiseasedBAI 
plot( years_HealthyBAI, difference_BAI, type="l", xlab="Anno", ylab="Differenza BAI (sani - malati)" )

# RWI (Ring Width Index)
## la crescita risponde diversamente a temperatura e precipitazioni?
### PC Gembloux
BeechTRWdetrend<-detrend(BeechTRW, method = "Spline", nyrs = 30)
BeechChron<-chron(BeechTRWdetrend,prefix = "AVG", biweight = TRUE, prewhiten = FALSE)
plot.crn(BeechChron)
### PC Sorsha
head(BeechTRW)
class(BeechTRW)
BeechTRW <- as.data.frame(BeechTRW)
rownames(BeechTRW) <- 1812:2018
BeechTRWdetrend <- detrend(BeechTRW, method = "Spline", nyrs = 30)
BeechChron<-chron(BeechTRWdetrend,prefix = "AVG", biweight = TRUE, prewhiten = FALSE)
plot.crn(BeechChron)
range(time(BeechChron))

plot(dcc(BeechChron, prec, selection = -6:9,method = "correlation",
           timespan = c(1930, 1990), var_names = "precipitation", boot = "std"))

plot(dcc(BeechChron, temp, selection = -6:9, method = "correlation",
           timespan = c(1930,2018), var_names = "temperature", boot = "std"))


### analisi clima-crescita
period <- c(1930, 1990)
month <- -6:9
plot(dcc(BeechChron,prec,selection = month, method = "correlation", 
            timespan = period, var_names = "precipitation", boot = "std"))

1) Periodo lungo (clima recente completo): 1930 - 1990
plot(dcc(BeechChron, prec, selection = -6:9, method = "correlation",
           timespan = c(1930,1990), var_names = "precipitation", boot = "std"))
plot(dcc(BeechChron, temp, selection = -6:9, method = "correlation",
         timespan = c(1930,1990), var_names = "temperature", boot = "std"))

2) Periodo storico più ampio: 1891 - 1990
plot(dcc(BeechChron, prec, selection = -6:9, method = "correlation",
           timespan = c(1891,1990), var_names = "precipitation", boot = "std"))
plot(dcc(BeechChron, temp, selection = -6:9, method = "correlation",
         timespan = c(1891,1990), var_names = "temperature", boot = "std"))

3) Periodo recente (cambiamento climatico): 1950 - 1990
plot(dcc(BeechChron, prec, selection = -6:9, method = "correlation",
           timespan = c(1950,1990), var_names = "precipitation", boot = "std"))
plot(dcc(BeechChron, temp, selection = -6:9, method = "correlation",
         timespan = c(1950,1990), var_names = "temperature", boot = "std"))


combiclim <-list(temp , prec)
ClimaBegin <- 1891
ClimaEnd <- 1990
plot(seascorr(BeechChron,combiclim,timespan=c(ClimaBegin,ClimaEnd),complete = 8, season_lengths = c(3),primary = 1,secondary = 2, ci = 0.05))
