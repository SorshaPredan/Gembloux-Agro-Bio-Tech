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


# BAI - BASAL AREA INCREMENT
## Do diseased trees grow less?
BeechBAI <- bai.in(TRW)

## TRW series 
matplot(
  as.numeric(rownames(TRW)),
  TRW,
  type = "l",
  xlab = "Year",
  ylab = "Tree-ring width"
)
## Standardization 
TRWdetrend <- detrend(
  TRW,
  method = "Spline",
  nyrs = 30
)
## Chronology CHRONOLOGY
BeechChron <- chron(
  TRWdetrend,
  prefix = "AVG",
  biweight = TRUE,
  prewhiten = FALSE
)
plot.crn(
  BeechChron
)
range(time(BeechChron))


# CLIMATE DATA
## PRECIPITATION
class(PrecSites)
head(PrecSites)
dim(PrecSites)

## Set column names
colnames(PrecSites) <- c(
  "year",
  "Jan", "Feb", "Mar", "Apr", "May", "Jun",
  "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
)
## Climate data are already numeric
PrecSites[, 2:13] <- lapply(
  PrecSites[, 2:13],
  as.numeric
)
## Check precipitation data
summary(PrecSites[, 2:13])
## Check for duplicated years
duplicati <- PrecSites$year[
  duplicated(PrecSites$year)
]
if (length(duplicati) > 0) {
  print(
    PrecSites[
      PrecSites$year %in% duplicati,
    ]
  )
}
## Remove duplicated years
PrecSites <- PrecSites[
  !duplicated(PrecSites$year),
]

# DENDROCLIMATIC CORRELATION - PRECIPITATION
plot(
  dcc(
    BeechChron,
    PrecSites,
    selection = -6:9,
    method = "correlation",
    timespan = c(1901, 2025),
    var_names = "precipitation",
    boot = "std"
  )
)

## TEMPERATURE
class(TempSites)
head(TempSites)
dim(TempSites)

## Set column names
colnames(TempSites) <- c(
  "year",
  "Jan", "Feb", "Mar", "Apr", "May", "Jun",
  "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
)
## Climate data are already numeric
TempSites[, 2:13] <- lapply(
  TempSites[, 2:13],
  as.numeric
)
## Check temperature data
summary(TempSites[, 2:13])
## Check for duplicated years
duplicati_temp <- TempSites$year[
  duplicated(TempSites$year)
]
if (length(duplicati_temp) > 0) {
  print(
    TempSites[
      TempSites$year %in% duplicati_temp,
    ]
  )
}
## Remove duplicated years
TempSites <- TempSites[
  !duplicated(TempSites$year),
]

# DENDROCLIMATIC CORRELATION - TEMPERATURE
plot(
  dcc(
    BeechChron,
    TempSites,
    selection = -6:9,
    method = "correlation",
    timespan = c(1901, 2025),
    var_names = "temperature",
    boot = "std"
  )
)



# CLIMATE-GROWTH CORRELATION
## Correlation TRW - Prec 
PrecCorr <- dcc(
  BeechChron,
  PrecSites,
  selection = -6:9,
  method = "correlation",
  timespan = c(1901, 2025),
  var_names = "precipitation",
  boot = "std"
)

## Correlation TRW - Temp
TempCorr <- dcc(
  BeechChron,
  TempSites,
  selection = -6:9,
  method = "correlation",
  timespan = c(1901, 2025),
  var_names = "temperature",
  boot = "std"
)



# COMBINED PEARSON CORRELATION PLOT
### Blue = precipitation
### Red  = temperature
### *    = significant correlation

## Extract correlation coefficients 
prec <- PrecCorr$coef$coef
temp <- TempCorr$coef$coef

## Define climate month labels 
mesi <- c(
  "Jun (-1)",
  "Jul (-1)",
  "Aug (-1)",
  "Sep (-1)",
  "Oct (-1)",
  "Nov (-1)",
  "Dec (-1)",
  "Jan",
  "Feb",
  "Mar",
  "Apr",
  "May",
  "Jun",
  "Jul",
  "Aug",
  "Sep"
)

## Extract significance
sig_prec <- PrecCorr$coef$significant
sig_temp <- TempCorr$coef$significant

## Create correlation matrix 
corr_matrix <- rbind(
  Precipitation = prec,
  Temperature = temp
)
colnames(corr_matrix) <- mesi


# Create Barplot 
bar_position <- barplot(
  corr_matrix,
  beside = TRUE,
  names.arg = mesi,
  col = c(
    "steelblue",
    "red"
  ),
  ylim = c(-0.75, 0.75),
  ylab = "Pearson r",
  xlab = "Climate month",
  las = 2,
  border = NA
)
abline(
  h = 0,
  lwd = 1.5
)
### Precipitation
text(
  x = bar_position[1, sig_prec],
  y = prec[sig_prec] +
    ifelse(
      prec[sig_prec] >= 0,
      0.04,
      -0.04
    ),
  labels = "*",
  cex = 1.5,
  col = "steelblue"
)
### Temperature
text(
  x = bar_position[2, sig_temp],
  y = temp[sig_temp] +
    ifelse(
      temp[sig_temp] >= 0,
      0.04,
      -0.04
    ),
  labels = "*",
  cex = 1.5,
  col = "red"
)
### Legend
legend(
  "topright",
  legend = c(
    "Precipitation",
    "Temperature"
  ),
  fill = c(
    "steelblue",
    "red"
  ),
  bty = "n"
)



# BEECH HEALTHY AND DISEASED #

### Do diseased trees show lower growth compared to healthy trees? ###
# Beech Malade/No Malade
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




# MOVING CLIMATE CORRELATIONS
## Healthy vs Diseased trees

# HEALTHY TREES - TEMPERATURE
Healthy_Temp_moving <- dcc(
  HealthyBeechChron,
  TempSites,
  selection = -6:9,
  method = "correlation",
  timespan = c(1901, 2025),
  var_names = "Temperature",
  boot = "std",
  dynamic = "moving"
)

plot(
  Healthy_Temp_moving,
  main = "Healthy trees - Moving correlation with Temperature"
)

# HEALTHY TREES - PRECIPITATION
Healthy_Prec_moving <- dcc(
  HealthyBeechChron,
  PrecSites,
  selection = -6:9,
  method = "correlation",
  timespan = c(1901, 2025),
  var_names = "Precipitation",
  boot = "std",
  dynamic = "moving"
)

plot(
  Healthy_Prec_moving,
  main = "Healthy trees - Moving correlation with Precipitation"
)



# DISEASED TREES - TEMPERATURE
Diseased_Temp_moving <- dcc(
  DiseasedBeechChron,
  TempSites,
  selection = -6:9,
  method = "correlation",
  timespan = c(1901, 2025),
  var_names = "Temperature",
  boot = "std",
  dynamic = "moving"
)

plot(
  Diseased_Temp_moving,
  main = "Diseased trees - Moving correlation with Temperature"
)

# DISEASED TREES - PRECIPITATION
Diseased_Prec_moving <- dcc(
  DiseasedBeechChron,
  PrecSites,
  selection = -6:9,
  method = "correlation",
  timespan = c(1901, 2025),
  var_names = "Precipitation",
  boot = "std",
  dynamic = "moving"
)

plot(
  Diseased_Prec_moving,
  main = "Diseased trees - Moving correlation with Precipitation"
)



# STATISTICAL COMPARISON
## Healthy vs Diseased climate response

### Prepare common years
Healthy_data <- data.frame(
  year = as.numeric(time(HealthyBeechChron)),
  Healthy = HealthyBeechChron$std
)
Diseased_data <- data.frame(
  year = as.numeric(time(DiseasedBeechChron)),
  Diseased = DiseasedBeechChron$std
)
Climate_data <- merge(
  Healthy_data,
  Diseased_data,
  by = "year"
)

## Function to compare Healthy and Diseased correlations
compare_correlations <- function(
    Healthy,
    Diseased,
    Climate,
    nboot = 1000
) {
  complete <- complete.cases(
    Healthy,
    Diseased,
    Climate
  )
  Healthy <- Healthy[complete]
  Diseased <- Diseased[complete]
  Climate <- Climate[complete]
  ### Observed correlations
  r_healthy <- cor(
    Healthy,
    Climate,
    method = "pearson"
  )
  r_diseased <- cor(
    Diseased,
    Climate,
    method = "pearson"
  )
  observed_difference <- r_healthy - r_diseased
  ### Bootstrap
  n <- length(Healthy)
  boot_difference <- numeric(nboot)
  set.seed(123)
  for (i in 1:nboot) {
    index <- sample(
      1:n,
      size = n,
      replace = TRUE
    )
    r_h <- cor(
      Healthy[index],
      Climate[index]
    )
    r_d <- cor(
      Diseased[index],
      Climate[index]
    )
    boot_difference[i] <- r_h - r_d
  }
  ### 95% confidence interval
  CI <- quantile(
    boot_difference,
    probs = c(0.025, 0.975),
    na.rm = TRUE
  )
  ### Two-sided p-value
  p_value <- 2 * min(
    mean(boot_difference <= 0),
    mean(boot_difference >= 0)
  )
  data.frame(
    Healthy_r = r_healthy,
    Diseased_r = r_diseased,
    Difference = observed_difference,
    CI_lower = CI[1],
    CI_upper = CI[2],
    p_value = p_value
  )
}


## TEMPERATURE
Temp_comparison <- data.frame()
for (month in c(
  "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
  "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL",
  "AUG", "SEP"
)) {
  ### Match month names between climate data and chronology
  if (month %in% names(TempSites)) {
    temp_month <- TempSites[
      TempSites$year %in% Climate_data$year,
      c("year", month)
    ]
    temp_month <- merge(
      Climate_data,
      temp_month,
      by = "year"
    )
    result <- compare_correlations(
      temp_month$Healthy,
      temp_month$Diseased,
      temp_month[[month]]
    )
    result$Climate <- "Temperature"
    result$Month <- month
    Temp_comparison <- rbind(
      Temp_comparison,
      result
    )
  }
}

## PRECIPITATION
Prec_comparison <- data.frame()
for (month in c(
  "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
  "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL",
  "AUG", "SEP"
)) {
  if (month %in% names(PrecSites)) {
    prec_month <- PrecSites[
      PrecSites$year %in% Climate_data$year,
      c("year", month)
    ]
    prec_month <- merge(
      Climate_data,
      prec_month,
      by = "year"
    )
    result <- compare_correlations(
      prec_month$Healthy,
      prec_month$Diseased,
      prec_month[[month]]
    )
    result$Climate <- "Precipitation"
    result$Month <- month
    Prec_comparison <- rbind(
      Prec_comparison,
      result
    )
  }
}


## COMBINE RESULTS
Climate_comparison <- rbind(
  Temp_comparison,
  Prec_comparison
)
### Reorder columns
Climate_comparison <- Climate_comparison[
  c(
    "Climate",
    "Month",
    "Healthy_r",
    "Diseased_r",
    "Difference",
    "CI_lower",
    "CI_upper",
    "p_value"
  )
]

### Round values
Climate_comparison$Healthy_r <-
  round(Climate_comparison$Healthy_r, 3)
Climate_comparison$Diseased_r <-
  round(Climate_comparison$Diseased_r, 3)
Climate_comparison$Difference <-
  round(Climate_comparison$Difference, 3)
Climate_comparison$CI_lower <-
  round(Climate_comparison$CI_lower, 3)
Climate_comparison$CI_upper <-
  round(Climate_comparison$CI_upper, 3)
Climate_comparison$p_value <-
  round(Climate_comparison$p_value, 3)


## Display results
Climate_comparison


# GRAPH: Healthy vs Diseased climate response
## Difference in correlation coefficients

### Order month
month_order <- c("Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
Climate_comparison$Month <- factor(
  Climate_comparison$Month,
  levels = month_order
)
### Order climate factor
Climate_comparison$Climate <- factor(
  Climate_comparison$Climate,
  levels = c("Temperature", "Precipitation")
)

## Graph
ggplot(
  Climate_comparison,
  aes(
    x = Month,
    y = Difference
  )
) +
  # linea dello zero
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black"
  ) +
  # intervalli di confidenza
  geom_errorbar(
    aes(
      ymin = CI_lower,
      ymax = CI_upper
    ),
    width = 0.15,
    color = "grey30"
  ) +
  # differenza
  geom_point(
    aes(
      color = p_value < 0.05
    ),
    size = 3
  ) +
  # separazione Temperature / Precipitation
  facet_wrap(
    ~ Climate,
    scales = "free_y"
  ) +
  # colori
  scale_color_manual(
    values = c(
      "FALSE" = "grey40",
      "TRUE" = "red"
    ),
    labels = c(
      "FALSE" = "Not significant",
      "TRUE" = "Significant"
    ),
    name = ""
  ) +
  labs(
    title = "Difference in climate-growth response",
    subtitle = "Healthy − Diseased",
    x = "Month",
    y = "Difference in Pearson correlation (r)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )



### BAI ###
HealthyBAI <- bai.in(HealthyTRW)
years_HealthyBAI <- as.numeric(row.names(HealthyBAI))
DiseasedBAI <- bai.in(DiseasedTRW)
years_DiseasedBAI <- as.numeric(row.names(DiseasedBAI))
## Annual mean chronology of healthy trees:
mean_HealthyBAI <- rowMeans(HealthyBAI, na.rm = TRUE)
## Annual mean chronology of diseased trees:
mean_DiseasedBAI <- rowMeans(DiseasedBAI, na.rm = TRUE)
plot(
  years_HealthyBAI,
  mean_HealthyBAI,
  type = "l",
  col = "blue",
  lwd = 2,
  xlab = "Year",
  ylab = expression("Basal Area Increment (cm"^2*")"),
  main = "Annual Basal Area Increment (BAI)"
)
lines(
  years_DiseasedBAI,
  mean_DiseasedBAI,
  col = "red",
  lwd = 2
)
legend(
  "topright",
  legend = c(
    "Healthy trees",
    "Diseased trees"
  ),
  col = c(
    "blue",
    "red"
  ),
  lwd = 2,
  bty = "n"
)

# Differences in growth between groups
difference_BAI <- mean_HealthyBAI - mean_DiseasedBAI 
plot( years_HealthyBAI, difference_BAI, type="l", xlab="Anno", ylab="Differenza BAI (sani - malati)" )                           
combiclim <-list(temp , prec)
ClimaBegin <- 1901
ClimaEnd <- 2025
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


## GRAPH: Temp
temp_coef <- result$coef[[1]]$primary$coef
temp_sig  <- result$coef[[1]]$primary$significant
temp_months <- c(
  "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
  "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul"
)
### Graph
bar_temp <- barplot(
  temp_coef,
  names.arg = temp_months,
  col = ifelse(temp_sig, "firebrick", "grey70"),
  ylim = c(-0.5, 0.5),
  ylab = "Pearson correlation (r)",
  xlab = "Season / month",
  main = "Beech growth – Temperature correlations",
  las = 2,
  border = NA
)
abline(h = 0, lwd = 1.5)
text(
  x = bar_temp[temp_sig],
  y = temp_coef[temp_sig] +
    0.04 * sign(temp_coef[temp_sig]),
  labels = "*",
  cex = 1.5
)
legend(
  "topright",
  legend = c("Significant", "Not significant"),
  fill = c("firebrick", "grey70"),
  bty = "n"
)


## GRAPH: Prec
prec_coef <- result$coef[[1]]$secondary$coef
prec_sig  <- result$coef[[1]]$secondary$significant
prec_months <- c(
  "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
  "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul"
)
bar_prec <- barplot(
  prec_coef,
  names.arg = prec_months,
  col = ifelse(prec_sig, "steelblue", "grey70"),
  ylim = c(-0.5, 0.5),
  ylab = "Pearson correlation (r)",
  xlab = "Season / month",
  main = "Beech growth – Precipitation correlations",
  las = 2,
  border = NA
)
abline(h = 0, lwd = 1.5)
if (any(prec_sig)) {
  text(
    x = bar_prec[prec_sig],
    y = prec_coef[prec_sig] +
      0.04 * sign(prec_coef[prec_sig]),
    labels = "*",
    cex = 1.5
  )
}
legend(
  "topright",
  legend = c("Significant", "Not significant"),
  fill = c("steelblue", "grey70"),
  bty = "n"
)
