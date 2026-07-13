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
BeechBAI<-bai.in(BeechTRW)
# rwi.stats(BirchTRW)
rwi.stats.running(BeechTRW)
corr.rwl.seg(BeechTRW)
## PC Gembloux
BeechTRWdetrend<-detrend(BeechTRW, method = "Spline", nyrs = 30)
BeechChron<-chron(BeechTRWdetrend,prefix = "AVG", biweight = TRUE, prewhiten = FALSE)
plot.crn(BeechChron)
## PC Sorsha
head(BeechTRW)
class(BeechTRW)
BeechTRW <- as.data.frame(BeechTRW)
BeechTRWdetrend <- detrend(BeechTRW, method = "Spline", nyrs = 30)
BeechChron<-chron(BeechTRWdetrend,prefix = "AVG", biweight = TRUE, prewhiten = FALSE)
plot.crn(BeechChron)
range(time(BeechChron))
range(as.numeric(rownames(BeechChron)))

plot(dcc(BeechChron, prec, selection = -6:9, method = "correlation",
           timespan = c(1930,1990), var_names = "precipitation", boot = "std"))

plot(dcc(BeechChron, temp, selection = -6:9, method = "correlation",
         timespan = c(1930,2018), var_names = "temperature", boot = "std"))

combiclim <-list(temp , prec)
plot(seascorr(BeechChron,combiclim,timespan=c(ClimaBegin,ClimaEnd),complete = 8, season_lengths = c(3),primary = 1,secondary = 2, ci = 0.05))

