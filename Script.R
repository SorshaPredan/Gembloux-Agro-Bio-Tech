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
TRW <- read.table(
   file.choose(),
   header = FALSE,
   fill = TRUE,
   stringsAsFactors = FALSE
 )
head(TRW)
dim(TRW)
campioni <- unique(TRW$V1)
serie <- list()
for(s in campioni){
   x <- TRW[TRW$V1 == s, ]
   anni <- c()
   valori <- c()
   for(i in 1:nrow(x)){
     anno_inizio <- x$V2[i]
     valori_riga <- as.numeric(x[i, 3:ncol(x)])
     valori_riga <- valori_riga[!is.na(valori_riga)]
     anni <- c(
       anni,
       anno_inizio:(anno_inizio + length(valori_riga) - 1)
     )
     valori <- c(valori, valori_riga)
   }
   serie[[s]] <- data.frame(
     anno = anni,
     valore = valori
   )
}
length(serie)
names(serie)[1:5]
head(serie[[1]])

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
# PC Gembloux
# NO PC Sorsha                     
TRW <- read.rwl("campioni.rwl")
head(TRW)
# Riprendere da qui                     
# Controllo qualità
rwi.stats(TRW)
rwi.stats.running(TRW)
corr.rwl.seg(TRW)

> library(openxlsx)
> write.xlsx(TRW, "C:/Users/user/Desktop/TIROCINIO FINALE/ANALISI DENDRO/pourSorsha/TRW.xlsx", rowNames = TRUE)
> ls()

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

plot(dcc(BeechChron, prec, selection = -6:9,method = "correlation",
           timespan = c(1930, 1990), var_names = "precipitation", boot = "std"))

plot(dcc(BeechChron, temp, selection = -6:9, method = "correlation",
           timespan = c(1930,2018), var_names = "temperature", boot = "std"))

# SITES
# Precipitation
class(PrecSites)
head(PrecSites)
dim(PrecSites)
colnames(PrecSites) <- c("year",
                         "Jan","Feb","Mar","Apr","May","Jun",
                         "Jul","Aug","Sep","Oct","Nov","Dec")
PrecSites <- PrecSites[, c("year",
                           "Jan","Feb","Mar","Apr","May","Jun",
                           "Jul","Aug","Sep","Oct","Nov","Dec")]
PrecSites$year <- 1901:2024
PrecSites <- PrecSites[, c("year", setdiff(names(PrecSites), "year"))]
PrecSites[,2:13] <- lapply(PrecSites[,2:13],
                            function(x) as.numeric(gsub("\\.", "", x)))
subset(PrecSites, year >= 1930 & year <= 1990) |> 
   summary()
plot(dcc(BeechChron, PrecSites, selection = -6:9, method = "correlation",
           timespan = c(1930,1990), var_names = "precipitation", boot = "std"))
# Temperature
class(TempSites)
head(TempSites)
dim(TempSites) 
colnames(TempSites) <- c("year",
                         "Jan","Feb","Mar","Apr","May","Jun",
                         "Jul","Aug","Sep","Oct","Nov","Dec")
TempSites[,2:13] <- lapply(TempSites[,2:13], function(x) {as.numeric(gsub("\\.", "", x))})
summary(TempSites[,2:13])  
TempSites[,2:13] <- TempSites[,2:13] / 1000000
any(is.na(TempSites))
diff(TempSites$year)
tail(TempSites$year)
TempSites <- TempSites[!duplicated(TempSites$year), ]                           
TempSites$year <- 1901:2023
plot(dcc(BeechChron, TempSites, selection = -6:9, method = "correlation",
           timespan = c(1930,1990), var_names = "temperature", boot = "std"))

### PLOT ###                       
Precipitation <- plot(dcc(BeechChron, PrecSites, selection = -6:9,method = "correlation",
           timespan = c(1930, 1990), var_names = "precipitation", boot = "std"))
Temperature <- plot(dcc(BeechChron, TempSites, selection = -6:9, method = "correlation",
           timespan = c(1930, 1990), var_names = "temperature", boot = "std"))                          

pdf("DCC_prec_temp.pdf", width = 8, height = 10)
par(mfrow = c(2,1))                          
plot(Precipitation, main="Precipitation 1930-1990")                          
plot(Temperature, main="Temperature 1930-1990")                           
dev.off()
                           
### analisi clima-crescita                           
month <- -6:9
# 1) Periodo lungo (serie climatica completa)
Prec1 <- plot(dcc(BeechChron, PrecSites, selection = month, method = "correlation",
           timespan = c(1902,2023), var_names = "precipitation", boot = "std"), main = "Precipitazioni 1901-2023")
Temp1 <- plot(dcc(BeechChron, TempSites, selection = month, method = "correlation",
           timespan = c(1902,2023), var_names = "temperature", boot = "std"), main = "Temperatura 1901-2023")
pdf("DCC_prec_temp.pdf", width = 8, height = 10)
par(mfrow = c(2,1))                          
plot(Prec1, main="Precipitation 1901-2023")                          
plot(Temp1, main="Temperature 1901-2023")                           
dev.off()                           
# 2) Periodo storico di confronto
Prec2 <- plot(dcc(BeechChron, PrecSites, selection = month, method = "correlation",
           timespan = c(1901,1990), var_names = "precipitation", boot = "std"), main = "Precipitazioni 1901-1990")
Temp2 <- plot(dcc(BeechChron, TempSites, selection = month, method = "correlation",
           timespan = c(1901,1990), var_names = "temperature", boot = "std"), main = "Temperatura 1901-1990")
pdf("DCC_prec_temp.pdf", width = 8, height = 10)
par(mfrow = c(2,1))                          
plot(Prec2, main="Precipitation 1901-1990")                          
plot(Temp2, main="Temperature 1901-1990")                           
dev.off()                           
# 3) Periodo recente (cambiamento climatico)
Prec3 <- plot(dcc(BeechChron, PrecSites, selection = month, method = "correlation",
           timespan = c(1950,2023), var_names = "precipitation", boot = "std"), main = "Precipitazioni 1950-2023")
Temp3 <- plot(dcc(BeechChron, TempSites, selection = month, method = "correlation",
           timespan = c(1950,2023), var_names = "temperature", boot = "std"), main = "Temperatura 1950-2023")
pdf("DCC_prec_temp.pdf", width = 8, height = 10)
par(mfrow = c(2,1))                          
plot(Prec3, main="Precipitation 1950-2023")                          
plot(Temp3, main="Temperature 1950-2023")                           
dev.off()

## Periodo evento 2009
Prec2009 <- plot(dcc(BeechChron, PrecSites, selection = 3:9, method = "correlation",
           timespan = c(2000,2010), var_names = "precipitation", boot = "std"))
Prec2009                           
Temp2009 <- plot(dcc(BeechChron, TempSites, selection = 3:9, method = "correlation", 
           timespan = c(2000,2010), var_names = "temperature", boot = "std")) 
Temp2009                           
pdf("DCC_prec_temp.pdf", width = 8, height = 10)
par(mfrow = c(2,1))                          
plot(Prec2009, main="Precipitation 2000-2010")                          
plot(Temp2009, main="Temperature 2000-2010")                           
dev.off()
## Periodo evento 2018-2019
Prec18 <- plot(dcc(BeechChron, PrecSites, selection = 3:9, method = "correlation",
           timespan = c(2010,2023), var_names = "precipitation", boot = "std"))
Temp18 <- plot(dcc(BeechChron, TempSites, selection = 3:9, method = "correlation",
          timespan = c(2010,2023), var_names = "temperature", boot = "std"))                   
pdf("DCC_prec_temp.pdf", width = 8, height = 10)
par(mfrow = c(2,1))                          
plot(Prec18, main="Precipitation 2010-2023")                          
plot(Temp18, main="Temperature 2010-2023")                           
dev.off()

# BAI Becch Malade/No Malade
## gli alberi malati presentano una crescita inferiore rispetto agli alberi sani?
HealthyTRW <- read_xlsx("TRWNoMalade.xlsx") 
HealthyTRW <- as.data.frame(HealthyTRW)
DiseasedTRW <- read_xlsx("TRWMalade.xlsx")
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
combiclim <-list(temp , prec)
ClimaBegin <- 1891
ClimaEnd <- 1990
plot(seascorr(BeechChron,combiclim,timespan=c(ClimaBegin,ClimaEnd),complete = 8, season_lengths = c(3),primary = 1,secondary = 2, ci = 0.05))
