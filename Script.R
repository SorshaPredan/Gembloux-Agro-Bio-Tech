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


prec<-read.table("Prec Chimay.txt", header = TRUE)
#temp <- read.table("tempBISHOP.txt", header = FALSE)
temp<-read.table("Temp Uccle.txt", header=TRUE)

#BirchTRW<- read_excel("TOTAL_TRW.xlsx")
BeechTRW<- read_xlsx("TRWBeech.xlsx")
#BeechTRW<-read_excel("TRWBeechMalade.xlsx")
BeechTRW[1]<-NULL
row.names(BeechTRW)<-1819:2025
rwi.stats(BeechTRW)

# BAI (Basal Area Increment = incrementO dell’area basale)
## gli alberi malati crescono meno?
BeechBAI<-bai.in(BeechTRW)     # Basal Area Increment = incrementdell’area basale
# rwi.stats(BirchTRW)
rwi.stats.running(BeechTRW)    # analizza come cambiano le statistiche della crescita nel tempo.
corr.rwl.seg(BeechTRW)         # controlla quanto gli alberi sono correlati tra loro

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
