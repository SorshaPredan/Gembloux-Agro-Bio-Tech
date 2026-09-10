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

# Climate Data
PrecSites<-read.table("Prec1901.txt", header = TRUE)
#temp <- read.table("tempBISHOP.txt", header = FALSE)
TempSites<-read.table("Temp1901.txt", header=TRUE)

# Preparation of the Tucson-format file for COFECHA analysis
## Import of the original dataset
TRW <- read_excel(
  file.choose(),
  col_names = FALSE
)
TRW <- as.data.frame(TRW)
head(TRW)
dim(TRW)
# File structure verification
## The first row contains the sample identifiers
## The first column contains the years (YEARS)
campioni <- as.character(TRW[1, -1])
anni <- TRW[-1, 1]

# Creation of the ring-width measurement matrix (RWL)
TRW <- TRW[-1, -1]
names(TRW) <- campioni
TRW[] <- lapply(TRW, as.numeric)
rownames(TRW) <- anni
TRW[TRW == 999] <- NA
class(TRW) <- c("rwl", "data.frame")
head(TRW)
dim(TRW)

# Export of the Tucson-format file for COFECHA analysis
write.rwl(
  TRW,
  "CAMPIONI.rwl",
  format = "tucson",
  long.names = TRUE
)

# DENDROCHRONOLOGICAL ANALYSIS
## Loading the cleaned Tucson-format file
### PC Gembloux
### NO PC Sorsha                     
TRW <- read.rwl("campioni.rwl")
head(TRW)
### Continue from this point                     
### Quality control
rwi.stats(TRW)
rwi.stats.running(TRW)
corr.rwl.seg(TRW)

library(openxlsx)
write.xlsx(TRW, "C:/Users/user/Desktop/TIROCINIO FINALE/ANALISI DENDRO/pourSorsha/TRW.xlsx", rowNames = TRUE)
ls()

# BAI (Basal Area Increment = incrementO dell’area basale)
## Do diseased trees grow less?
BeechBAI <- bai.in(TRW)
### Plot of raw tree-ring series
matplot(
  as.numeric(rownames(TRW)),
  TRW,
  type="l",
  xlab="Year",
  ylab="Tree ring width"
)                      
### Standardization
TRWdetrend<-detrend(TRW, method = "Spline", nyrs = 30)
### Chronology                      
BeechChron<-chron(TRWdetrend,prefix = "AVG", biweight = TRUE, prewhiten = FALSE)
plot.crn(BeechChron)                      
range(time(BeechChron))

# SITES
## Precipitation
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

## Temperature
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

# CLIMATE-GROWTH CORRELATION
## Correlation TRW-Precipitation 
PrecCorr <- dcc(BeechChron, PrecSites, selection = -6:9, method = "correlation",
                  timespan = c(1930,1990), var_names = "precipitation",boot = "std")
## Correlation TRW-Temperature
TempCorr <- dcc(BeechChron,TempSites,selection = -6:9,method = "correlation",
                  timespan = c(1930,1990), var_names = "temperature", boot = "std")

# Combined Pearson correlation plot
### Blue = precipitation
### Red = temperature
## Extraction of Pearson correlation coefficients
prec <- PrecCorr$coef$coef
temp <- TempCorr$coef$coef
### month
mesi <- PrecCorr$coef$month
### matrix preparation for the plot
corr_matrix <- rbind(
  prec,
  temp
)
### plot
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
### line at zero
abline(h = 0, lwd = 2)
### legend
legend(
  "topright",
  legend = c("Precipitation", "Temperature"),
  fill = c("steelblue", "red"),
  bty = "n"
)
# Statistical significance of correlations
## * = significant correlations (DCC bootstrap)
### extraction of significance values
sig_prec <- PrecCorr$coef$significant
sig_temp <- TempCorr$coef$significant
### creation of the barplot
bar_position <- barplot(
  corr_matrix,
  beside = TRUE,
  names.arg = colnames(corr_matrix),
  col = c("steelblue", "red"),
  ylim = c(-0.75, 0.75),
  las = 2
)
### extraction of the positions of the two bar groups
pos_prec <- bar_position[1, ]
pos_temp <- bar_position[2, ]
### addition of significance markers for precipitation
text(
  pos_prec[sig_prec],
  prec[sig_prec] + 0.05 * sign(prec[sig_prec]),
  "*",
  cex = 1.5,
  col = "steelblue"
)
### addition of significance markers for temperature
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
                           
### Climate–growth relationship analysis ###                       
month <- -6:9
# 1) Long-term period (complete climate series)
Prec1 <- plot(dcc(BeechChron, PrecSites, selection = month, method = "correlation",
           timespan = c(1902,2023), var_names = "precipitation", boot = "std"), main = "Precipitazioni 1901-2023")
Temp1 <- plot(dcc(BeechChron, TempSites, selection = month, method = "correlation",
           timespan = c(1902,2023), var_names = "temperature", boot = "std"), main = "Temperatura 1901-2023")
pdf("DCC_prec_temp.pdf", width = 8, height = 10)
par(mfrow = c(2,1))                          
plot(Prec1, main="Precipitation 1901-2023")                          
plot(Temp1, main="Temperature 1901-2023")                           
dev.off()                           
# 2) Historical reference period
Prec2 <- plot(dcc(BeechChron, PrecSites, selection = month, method = "correlation",
           timespan = c(1901,1990), var_names = "precipitation", boot = "std"), main = "Precipitazioni 1901-1990")
Temp2 <- plot(dcc(BeechChron, TempSites, selection = month, method = "correlation",
           timespan = c(1901,1990), var_names = "temperature", boot = "std"), main = "Temperatura 1901-1990")
pdf("DCC_prec_temp.pdf", width = 8, height = 10)
par(mfrow = c(2,1))                          
plot(Prec2, main="Precipitation 1901-1990")                          
plot(Temp2, main="Temperature 1901-1990")                           
dev.off()                           
# 3) Recent period under climate change
Prec3 <- plot(dcc(BeechChron, PrecSites, selection = month, method = "correlation",
           timespan = c(1950,2023), var_names = "precipitation", boot = "std"), main = "Precipitazioni 1950-2023")
Temp3 <- plot(dcc(BeechChron, TempSites, selection = month, method = "correlation",
           timespan = c(1950,2023), var_names = "temperature", boot = "std"), main = "Temperatura 1950-2023")
pdf("DCC_prec_temp.pdf", width = 8, height = 10)
par(mfrow = c(2,1))                          
plot(Prec3, main="Precipitation 1950-2023")                          
plot(Temp3, main="Temperature 1950-2023")                           
dev.off()


### Do diseased trees show lower growth compared to healthy trees? ###
# Becch Malade/No Malade
HealthyTRW <- read_xlsx("TRWNoMalade.xlsx") 
HealthyTRW <- as.data.frame(HealthyTRW)
DiseasedTRW <- read_xlsx("TRWMalade.xlsx")
DiseasedTRW <- as.data.frame(DiseasedTRW)
HealthyTRW[1] <- NULL
DiseasedTRW[1] <- NULL
row.names(HealthyTRW)<-1800:2025
row.names(DiseasedTRW)<-1800:2025

## Tree-ring quality control
rwi.stats(HealthyTRW)
rwi.stats.running(HealthyTRW)
corr.rwl.seg(HealthyTRW)
## Tree-ring quality control
rwi.stats(DiseasedTRW)
rwi.stats.running(DiseasedTRW)
corr.rwl.seg(DiseasedTRW) 
## Tree-ring standardization
HealthyTRWdetrend<-detrend(HealthyTRW, method = "Spline", nyrs = 30)
## Tree-ring chronology                   
HealthyBeechChron<-chron(HealthyTRWdetrend,prefix = "AVG", biweight = TRUE, prewhiten = FALSE)
plot.crn(HealthyBeechChron)                      
range(time(HealthyBeechChron))   
## Tree-ring standardization                          
DiseasedTRWdetrend<-detrend(DiseasedTRW, method = "Spline", nyrs = 30)
## Tree-ring chronology                    
DiseasedBeechChron<-chron(DiseasedTRWdetrend,prefix = "AVG", biweight = TRUE, prewhiten = FALSE)
plot.crn(DiseasedBeechChron)                      
range(time(DiseasedBeechChron))

# Creation of an empty plot using the common time range
plot(time(HealthyBeechChron), HealthyBeechChron$std,
     type = "l",
     col = "blue",
     lwd = 2,
     ylim = range(c(HealthyBeechChron$std, DiseasedBeechChron$std), na.rm = TRUE),
     xlab = "Year",
     ylab = "Ring width index",
     main = "Beech chronology: Healthy vs Diseased")
### addition of the diseased tree-ring chronology
lines(time(DiseasedBeechChron), DiseasedBeechChron$std,
      col = "red",
      lwd = 2)
### legend
legend("topright",
       legend = c("Healthy", "Diseased"),
       col = c("blue", "red"),
       lwd = 2) 

### BAI ###
HealthyBAI <- bai.in(HealthyTRW)
years_HealthyBAI <- as.numeric(row.names(HealthyBAI))
DiseasedBAI <- bai.in(DiseasedTRW)
years_DiseasedBAI <- as.numeric(row.names(DiseasedBAI))
## Annual mean chronology of healthy trees:
mean_HealthyBAI <- rowMeans(HealthyBAI, na.rm = TRUE)
## Annual mean chronology of diseased trees:
mean_DiseasedBAI <- rowMeans(DiseasedBAI, na.rm = TRUE)
plot( years_HealthyBAI, mean_HealthyBAI, type="l", col="blue", lwd=2, xlab="Year", ylab="BAI (cm²/anno)" ) 
      lines( years_DiseasedBAI, mean_DiseasedBAI, col="red", lwd=2 ) 
     legend( "topright", legend=c("Beech No Malade", "Beech Malade"), col=c("blue","red"), lwd=2 )

# Differences in growth between groups
difference_BAI <- mean_HealthyBAI - mean_DiseasedBAI 
plot( years_HealthyBAI, difference_BAI, type="l", xlab="Anno", ylab="Differenza BAI (sani - malati)" )                           
combiclim <-list(temp , prec)
ClimaBegin <- 1902
ClimaEnd <- 2023
### Creation of the combined climate dataset
combiclim <- list(
  PrecSites,
  TempSites
)
### Seasonal correlation analysis
result <- seascorr(
  BeechChron,
  combiclim,
  timespan = c(ClimaBegin, ClimaEnd),
  complete = 8,
  season_lengths = c(5),
  primary = 1,
  secondary = 2,
  ci = 0.05
)
### Plot customization
plot(result) +
  scale_fill_manual(
    values = c("grey85", "firebrick")
  ) +
  ggtitle("Seasonal climate-growth correlations") +
  theme_bw() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 14,
      face = "bold"
    ),
    legend.position = "bottom",
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )
