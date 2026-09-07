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

### SPEI SUMMER ###
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
## SPEI classification:
### 0 = normal conditions
### negative values indicate drought
### -1 = moderate drought
### -1.5 = severe drought
### -2 = extreme drought
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

### SEA (Superposed Epoch Analysis) ###
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

# PC Sorsha
XCT_folder <- "C:/Users/user/Desktop/TIROCINIO FINALE/DENSITA/XCT/CORES"

# EWD DENSITY
EWD_sea_data <- data.frame(
  EWD = EWDChron$std
)
rownames(EWD_sea_data) <- rownames(EWDChron)
class(EWD_sea_data) <- c("rwl","data.frame")
EWD_SEA <- sea(
  EWD_sea_data,
  key = drought_events,
  lag = 3,
  resample = 1000
)
plot(EWD_SEA)

plot(
  EWD_SEA$lag,
  EWD_SEA$se,
  type = "b",
  pch = 16,
  col = "darkgreen",
  lwd = 2,
  cex = 1.2,
  ylim = c(-0.8, -0.5),
  main = "SEA - EWD response to extreme drought",
  xlab = "Years relative to drought event",
  ylab = "SEA response"
)
abline(h = 0, lty = 2, col = "grey40")
abline(v = 0, lty = 2, col = "red")

# LWD DENSITY
LWD_sea_data <- data.frame(
  LWD = LWDChron$std
)
rownames(LWD_sea_data) <- rownames(LWDChron)
class(LWD_sea_data) <- c("rwl","data.frame")
LWD_SEA <- sea(
  LWD_sea_data,
  key = drought_events,
  lag = 3,
  resample = 1000
)
plot(LWD_SEA)

plot(
  LWD_SEA$lag,
  LWD_SEA$se,
  type = "b",
  pch = 16,
  col = "navy",
  lwd = 2,
  cex = 1.2,
  main = "SEA - LWD response to extreme drought",
  xlab = "Years relative to drought event",
  ylab = "SEA response"
)
abline(h = 0, lty = 2, col = "grey40")
abline(v = 0, lty = 2, col = "red")

# MXD DENSITY
MXD_sea_data <- data.frame(
  MXD = MXDChron$std
)
rownames(MXD_sea_data) <- rownames(MXDChron)
class(MXD_sea_data) <- c("rwl","data.frame")
MXD_SEA <- sea(
  MXD_sea_data,
  key = drought_events,
  lag = 3,
  resample = 1000
)
plot(MXD_SEA)

plot(
  MXD_SEA$lag,
  MXD_SEA$se,
  type = "b",
  pch = 16,
  col = "steelblue",
  lwd = 2,
  cex = 1.2,
  main = "SEA - MXD response to extreme drought",
  xlab = "Years relative to drought event",
  ylab = "SEA response",
  ylim = range(MXD_SEA$se) + c(-0.05, 0.05)
)
abline(h = 0, lty = 2, col = "grey40")
abline(v = 0, lty = 2, col = "red")

## Superposed Epoch Analysis of EWD, LWD, and MXD
SEA_all <- bind_rows(
  data.frame(
    lag = EWD_SEA$lag,
    response = EWD_SEA$se,
    p = EWD_SEA$p,
    parameter = "EWD"
  ),
  data.frame(
    lag = LWD_SEA$lag,
    response = LWD_SEA$se,
    p = LWD_SEA$p,
    parameter = "LWD"
  ),
  data.frame(
    lag = MXD_SEA$lag,
    response = MXD_SEA$se,
    p = MXD_SEA$p,
    parameter = "MXD"
  )
)
SEA_all <- SEA_all %>%
  mutate(
    sig = case_when(
      p < 0.01 ~ "**",
      p < 0.05 ~ "*",
      TRUE ~ ""
    )
  )
head(SEA_all)

### plot
library(ggplot2)
p_SEA <- ggplot(
  SEA_all,
  aes(
    x = lag,
    y = response,
    colour = parameter
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey40"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "red"
  ) +
  geom_text(
    aes(
      label = sig,
      y = response + 0.05
    ),
    colour = "black",
    size = 5
  ) +
  scale_colour_manual(
    values = c(
      EWD = "darkgreen",
      LWD = "orange",
      MXD = "red"
    )
  ) +
 labs(
  x = "Years relative to drought event",
  y = "Mean standardized response",
  colour = "Parameter",
  title = "Superposed Epoch Analysis – Extreme drought events",
  subtitle = "Drought years: 1922, 1927, 1953, 1965, 1972, 1980, 1987, 2007, 2012",
  caption = "* p < 0.05; ** p < 0.01"
) +
  theme_classic() +
theme(
  plot.caption = element_text(
    hjust = 0,
    size = 10
  )
)

p_SEA
