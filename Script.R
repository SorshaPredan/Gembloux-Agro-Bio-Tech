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

# GRAFICO COMBINATO PEARSON r
### Blu = precipitazioni
### Rosso = temperatura
# estrazione coefficienti Pearson
prec <- PrecCorr$coef$coef
temp <- TempCorr$coef$coef
# mesi
mesi <- PrecCorr$coef$month
# matrice per il grafico
corr_matrix <- rbind(
  prec,
  temp
)
# grafico
barplot(
  rbind(prec, temp),
  beside = TRUE,
  names.arg = mesi,
  col = c("steelblue", "red"),
  ylim = c(-0.75,0.75),
  ylab = "Pearson r",
  xlab = "Month",
  las = 2
)
# linea dello zero
abline(h = 0, lwd = 2)
# legenda
legend(
  "topright",
  legend = c("Precipitation", "Temperature"),
  fill = c("steelblue", "red"),
  bty = "n"
)
# Significatività delle correlazioni
## * = correlazione significativa (bootstrap dcc)
# recupero significatività
sig_prec <- PrecCorr$coef$significant
sig_temp <- TempCorr$coef$significant
# recupero posizioni delle barre
bar_position <- barplot(
  corr_matrix,
  beside = TRUE,
  plot = FALSE
)
pos_prec <- bar_position[1,]
pos_temp <- bar_position[2,]
# aggiunta degli asterischi
text(
  pos_prec[sig_prec],
  prec[sig_prec] + 0.05 * sign(prec[sig_prec]),
  "*",
  cex = 1.5,
  col = "steelblue"
)
text(
  pos_temp[sig_temp],
  temp[sig_temp] + 0.05 * sign(temp[sig_temp]),
  "*",
  cex = 1.5,
  col = "red"
)

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

# Becch Malade/No Malade
## gli alberi malati presentano una crescita inferiore rispetto agli alberi sani?
HealthyTRW <- read_xlsx("TRWNoMalade.xlsx") 
HealthyTRW <- as.data.frame(HealthyTRW)
DiseasedTRW <- read_xlsx("TRWMalade.xlsx")
DiseasedTRW <- as.data.frame(DiseasedTRW)
HealthyTRW[1] <- NULL
DiseasedTRW[1] <- NULL
row.names(HealthyTRW)<-1800:2025
row.names(DiseasedTRW)<-1800:2025

# Controllo qualità
rwi.stats(HealthyTRW)
rwi.stats.running(HealthyTRW)
corr.rwl.seg(HealthyTRW)
# Controllo qualità
rwi.stats(DiseasedTRW)
rwi.stats.running(DiseasedTRW)
corr.rwl.seg(DiseasedTRW) 
# Standardizzazione
HealthyTRWdetrend<-detrend(HealthyTRW, method = "Spline", nyrs = 30)
# Cronologia                      
HealthyBeechChron<-chron(HealthyTRWdetrend,prefix = "AVG", biweight = TRUE, prewhiten = FALSE)
plot.crn(HealthyBeechChron)                      
range(time(HealthyBeechChron))   
# Standardizzazione                           
DiseasedTRWdetrend<-detrend(DiseasedTRW, method = "Spline", nyrs = 30)
# Cronologia                      
DiseasedBeechChron<-chron(DiseasedTRWdetrend,prefix = "AVG", biweight = TRUE, prewhiten = FALSE)
plot.crn(DiseasedBeechChron)                      
range(time(DiseasedBeechChron))

# Creazione del grafico vuoto usando il range temporale comune
plot(time(HealthyBeechChron), HealthyBeechChron$std,
     type = "l",
     col = "blue",
     lwd = 2,
     ylim = range(c(HealthyBeechChron$std, DiseasedBeechChron$std), na.rm = TRUE),
     xlab = "Year",
     ylab = "Ring width index",
     main = "Beech chronology: Healthy vs Diseased")
# Aggiunta della cronologia Diseased
lines(time(DiseasedBeechChron), DiseasedBeechChron$std,
      col = "red",
      lwd = 2)
# Legenda
legend("topright",
       legend = c("Healthy", "Diseased"),
       col = c("blue", "red"),
       lwd = 2) 

### BAI ###
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
