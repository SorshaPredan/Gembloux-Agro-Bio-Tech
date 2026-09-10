# CLIMATE - RWI ANALYSIS
# ANNUAL MOVING CORRELATION + DCC
## Period: 1901-2025


library(readxl)
library(dplR)
library(openxlsx)
library(ggplot2)
library(dplyr)
library(tidyr)

### DCC
## install.packages("rugarch")
## install.packages("rmgarch")
library(rugarch)
library(rmgarch)


## Climate data
PrecSites <- read.table(
  "Prec1901.txt",
  header = TRUE
)
TempSites <- read.table(
  "Temp1901.txt",
  header = TRUE
)
head(PrecSites)
head(TempSites)
str(PrecSites)
str(TempSites)

## Import data
TRW <- read_excel(
  file.choose(),
  col_names = FALSE
)
TRW <- as.data.frame(TRW)
head(TRW)
dim(TRW)

## Preparation of RWL Matrix 4. PREPARATION OF RWL MATRIX
### First row = sample identifiers
campioni <- as.character(TRW[1, -1])
### First column = years
anni <- TRW[-1, 1]
### Remove first row and first column
TRW <- TRW[-1, -1]
### Assign sample names
names(TRW) <- campioni
### Convert to numeric
TRW[] <- lapply(
  TRW,
  as.numeric
)
### Assign years
rownames(TRW) <- anni
### 999 = missing value
TRW[TRW == 999] <- NA
### Set dplR class
class(TRW) <- c(
  "rwl",
  "data.frame"
)
head(TRW)
dim(TRW)

## Quality control
rwi.stats(TRW)
rwi.stats.running(TRW)
corr.rwl.seg(TRW)

## Save RWL 6. SAVE RWL
write.xlsx(
  TRW,
  "TRW.xlsx",
  rowNames = TRUE
)

## Detrending
### Spline detrending
RWI <- detrend(
  TRW,
  method = "Spline"
)
head(RWI)
dim(RWI)

## Build RWI Chronology 8. BUILD RWI CHRONOLOGY
### Robust chronology from the individual RWI series
RWI_chron <- chron(
  RWI
)
head(RWI_chron)
dim(RWI_chron)
### Inspect column names
colnames(RWI_chron)



## Prepare annual RWI series 
### Convert row names into a numeric year column
RWI_chron$year <- as.numeric(
  rownames(RWI_chron)
)
### Rename the chronology column
### In dplR this is normally "std"
### We explicitly identify it to avoid relying on column position.
rwi_column <- setdiff(
  colnames(RWI_chron),
  c("year", "sample.depth", "res", "res.se", "std.err")
)
### If "std" exists, use it
if ("std" %in% colnames(RWI_chron)) {
  
  RWI_annual <- RWI_chron %>%
    select(
      year,
      RWI = std
    )  
} else {
  
### Fallback: first numeric chronology column
  RWI_annual <- RWI_chron %>%
    select(
      year,
      RWI = all_of(rwi_column[1])
    )
}


## Keep common period 1901-2025
RWI_annual <- RWI_annual %>%
  filter(
    year >= 1901,
    year <= 2025
  )


## Prepare annual climate data 
# TEMPERATURE
Temperature_annual <- data.frame(
  year = TempSites$Year,
  temperature = rowMeans(
    TempSites[, paste0("X", 1:12)],
    na.rm = TRUE
  )
)

# PRECIPITATION
Precipitation_annual <- data.frame(
  year = PrecSites$Year,
  precipitation = rowSums(
    PrecSites[, paste0("X", 1:12)],
    na.rm = TRUE
  )
)
head(Temperature_annual)
head(Precipitation_annual)
str(Temperature_annual)
str(Precipitation_annual)

## Marge RWI + Temp + Prec 
climate_rwi <- merge(
  RWI_annual,
  Temperature_annual,
  by = "year",
  all = FALSE
)
climate_rwi <- merge(
  climate_rwi,
  Precipitation_annual,
  by = "year",
  all = FALSE
)
### Sort by year
climate_rwi <- climate_rwi[
  order(climate_rwi$year),
]
### Keep only the common period
climate_rwi <- climate_rwi[
  climate_rwi$year >= 1901 &
  climate_rwi$year <= 2025,
]
head(climate_rwi)
tail(climate_rwi)
dim(climate_rwi)
summary(climate_rwi)
colSums(is.na(climate_rwi))

## 25-year moving correlation 
moving_cor <- function(x, y, window = 25) {
  result <- rep(
    NA_real_,
    length(x)
  )
  for (i in seq_along(x)) {
    if (i < window) {
      next
    }
    index <- (i - window + 1):i
    result[i] <- cor(
      x[index],
      y[index],
      method = "pearson"
    )
  }
  return(result)
}

## RWI - TEMPERATURE
climate_rwi$cor_temp_25yr <- moving_cor(
  climate_rwi$RWI,
  climate_rwi$temperature,
  window = 25
)

## RWI - PRECIPITATION
climate_rwi$cor_prec_25yr <- moving_cor(
  climate_rwi$RWI,
  climate_rwi$precipitation,
  window = 25
)
head(climate_rwi,30)
tail(climate_rwi)

## Data for heatmap 14. PREPARE DATA FOR HEATMAP
grafico <- climate_rwi[
  ,
  c(
    "year",
    "cor_temp_25yr",
    "cor_prec_25yr"
  )
]
names(grafico) <- c(
  "year",
  "Temperature",
  "Precipitation"
)
### Convert from wide to long format
grafico <- reshape(
  grafico,
  varying = c(
    "Temperature",
    "Precipitation"
  ),
  v.names = "correlation",
  timevar = "variable",
  times = c(
    "Temperature",
    "Precipitation"
  ),
  direction = "long"
)
rownames(grafico) <- NULL
### Order of variables
grafico$variable <- factor(
  grafico$variable,
  levels = c(
    "Precipitation",
    "Temperature"
  )
)

## Annual 25-year moving correlation heatmap
# CORRELATION CLASSES
grafico$classe <- cut(
  grafico$correlation,
  breaks = c(
    -1,
    -0.75,
    -0.50,
    -0.25,
    0,
    0.25,
    0.50,
    0.75,
    1
  ),
  labels = c(
    "-1.00 – -0.75",
    "-0.75 – -0.50",
    "-0.50 – -0.25",
    "-0.25 – 0.00",
    "0.00 – 0.25",
    "0.25 – 0.50",
    "0.50 – 0.75",
    "0.75 – 1.00"
  ),
  include.lowest = TRUE
)


# HEATMAP
p <- ggplot(
  grafico,
  aes(
    x = year,
    y = variable,
    fill = classe
  )
) +
  geom_tile(
    width = 1,
    height = 0.95,
    colour = "black",
    linewidth = 0.45
  ) +
  
  ### high-contrast colours
  scale_fill_manual(
    values = c(
      "-1.00 – -0.75" = "#001F8B",  # blu molto scuro
      "-0.75 – -0.50" = "#0057D9",  # blu intenso
      "-0.50 – -0.25" = "#00A6D6",  # azzurro/ciano
      "-0.25 – 0.00"  = "#5CE1E6",  # ciano chiaro
      
      "0.00 – 0.25"   = "#FFF200",  # giallo acceso
      "0.25 – 0.50"   = "#FFB000",  # giallo/arancio
      "0.50 – 0.75"   = "#FF5A00",  # arancio intenso
      "0.75 – 1.00"   = "#D9003F"   # rosso/magenta
    ),
    na.value = "#BDBDBD",
    drop = FALSE,
    name = "Pearson r"
  ) +

  ### X axis
  scale_x_continuous(
    breaks = seq(
      1905,
      2025,
      by = 5
    ),
    expand = c(
      0,
      0
    )
  ) +

  ### Y axis
  scale_y_discrete(
    labels = c(
      "Precipitation" = "Precipitation",
      "Temperature" = "Temperature"
    ),
    expand = c(
      0,
      0
    )
  ) +
  
  ### titles TITLES
  labs(
    title =
      "25-Year Moving Correlation between RWI and Climate",
    subtitle =
      "Annual temperature and precipitation | 1901–2025",
    x =
      "Year",
    y =
      NULL
  ) +
  
  ### theme
  theme_minimal() +
  theme(
    panel.grid =
      element_blank(),
    axis.text.y =
      element_text(
        size = 11,
        colour = "black",
        face = "bold"
      ),
    axis.text.x =
      element_text(
        size = 8,
        colour = "black",
        angle = 45,
        hjust = 1
      ),
    axis.title.x =
      element_text(
        size = 10,
        face = "bold"
      ),
    plot.title =
      element_text(
        size = 14,
        face = "bold",
        hjust = 0.5
      ),
    plot.subtitle =
      element_text(
        size = 10,
        hjust = 0.5
      ),
    legend.title =
      element_text(
        size = 9,
        face = "bold"
      ),
    legend.text =
      element_text(
        size = 7.5
      ),
    legend.key.height =
      unit(
        0.65,
        "cm"
      ),
    plot.margin =
      margin(
        10,
        15,
        10,
        10
      )
  )

# Display graph
print(p)



# Moving Correlation
## TEMPERATURE-RWI vs PRECIPITATION-RWI

library(ggplot2)

grafico_linee <- rbind(
  data.frame(
    year = climate_rwi$year,
    correlation = climate_rwi$cor_temp_25yr,
    variable = "Temperature"
  ),
  data.frame(
    year = climate_rwi$year,
    correlation = climate_rwi$cor_prec_25yr,
    variable = "Precipitation"
  )
)

### order of panel 
grafico_linee$variable <- factor(
  grafico_linee$variable,
  levels = c(
    "Temperature",
    "Precipitation"
  )
)



# GRAPH
p_linee <- ggplot(
  grafico_linee,
  aes(
    x = year,
    y = correlation
  )
) +
  geom_hline(
    yintercept = 0,
    colour = "black",
    linewidth = 0.5
  ) +
  geom_line(
    colour = "#D9003F",
    linewidth = 1
  ) +
  facet_wrap(
    ~ variable,
    ncol = 1,
    scales = "fixed"
  ) +
  scale_x_continuous(
    breaks = seq(
      1925,
      2025,
      by = 5
    ),
    expand = c(
      0,
      0
    )
  ) +
  scale_y_continuous(
    limits = c(
      -1,
      1
    ),
    breaks = seq(
      -1,
      1,
      by = 0.25
    )
  ) +
  labs(
    title =
      "25-Year Moving Correlation between RWI and Climate",
    subtitle =
      "Annual temperature and precipitation | 1901–2025",
    x =
      "Year",
    y =
      "Pearson r"
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor =
      element_blank(),
    panel.grid.major.x =
      element_blank(),
    strip.text =
      element_text(
        size = 11,
        face = "bold",
        colour = "black"
      ),
    strip.background =
      element_rect(
        fill = "#F2F2F2",
        colour = "black",
        linewidth = 0.4
      ),
    axis.text.x =
      element_text(
        size = 8,
        colour = "black",
        angle = 45,
        hjust = 1
      ),
    axis.text.y =
      element_text(
        size = 8,
        colour = "black"
      ),
    axis.title.x =
      element_text(
        size = 10,
        face = "bold"
      ),
    axis.title.y =
      element_text(
        size = 10,
        face = "bold"
      ),
    plot.title =
      element_text(
        size = 14,
        face = "bold",
        hjust = 0.5
      ),
    plot.subtitle =
      element_text(
        size = 10,
        hjust = 0.5
      ),
    plot.margin =
      margin(
        10,
        15,
        10,
        10
      )
  )


# Display
print(p_linee)



# DCC-GARCH ANALYSIS 

library(rmgarch)
library(rugarch)


head(climate_rwi)
summary(climate_rwi[, c(
  "RWI",
  "temperature",
  "precipitation"
)])
colSums(is.na(climate_rwi[, c(
  "RWI",
  "temperature",
  "precipitation"
)]))

### Creation of the series 
RWI_dcc <- climate_rwi$RWI
Temperature_dcc <- climate_rwi$temperature
Precipitation_dcc <- climate_rwi$precipitation

### Conversion to time series 
RWI_dcc <- ts(
  RWI_dcc,
  start = 1901,
  frequency = 1
)
Temperature_dcc <- ts(
  Temperature_dcc,
  start = 1901,
  frequency = 1
)
Precipitation_dcc <- ts(
  Precipitation_dcc,
  start = 1901,
  frequency = 1
)

start(RWI_dcc)
end(RWI_dcc)
start(Temperature_dcc)
end(Temperature_dcc)
start(Precipitation_dcc)
end(Precipitation_dcc)
length(RWI_dcc)
length(Temperature_dcc)
length(Precipitation_dcc)


# DCC-GARCH: RWI - TEMPERATURE
data_dcc_temp <- cbind(
  RWI = as.numeric(RWI_dcc),
  Temperature = as.numeric(Temperature_dcc)
)
head(data_dcc_temp)
dim(data_dcc_temp)
colSums(is.na(data_dcc_temp))

## Univariate GARCH SPECIFICATION
spec_garch <- ugarchspec(
  variance.model = list(
    model = "sGARCH",
    garchOrder = c(1, 1)
  ),
  mean.model = list(
    armaOrder = c(1, 0),
    include.mean = TRUE
  ),
  distribution.model = "norm"
)

### DCC specification
spec_dcc_temp <- dccspec(
  
  uspec = multispec(
    replicate(
      2,
      spec_garch
    )
  ),
  dccOrder = c(1, 1),
  distribution = "mvnorm"
)


### Fit DCC model 20.4 FIT DCC MODEL
dcc_temp <- dccfit(
  spec_dcc_temp,
  data = data_dcc_temp,
  fit.control = list(
    eval.se = TRUE
  )
)
show(dcc_temp)


# DCC-GARCH: RWI - PRECIPITATION
data_dcc_prec <- cbind(
  RWI = as.numeric(RWI_dcc),
  Precipitation = as.numeric(Precipitation_dcc)
)
head(data_dcc_prec)
dim(data_dcc_prec)
colSums(is.na(data_dcc_prec))

## DCC specification 21.2 DCC SPECIFICATION
spec_dcc_prec <- dccspec(
  uspec = multispec(
    replicate(
      2,
      spec_garch
    )
  ),
  dccOrder = c(1, 1),
  distribution = "mvnorm"
)

### Fit DCC model 21.3 FIT DCC MODEL
dcc_prec <- dccfit( 
  spec_dcc_prec,
  data = data_dcc_prec,
  fit.control = list(
    eval.se = TRUE
  )
)
show(dcc_prec)


## Extract dynamic conditional correlation 
cor_dcc_temp <- rcor(dcc_temp)
cor_dcc_temp <- as.numeric(
  cor_dcc_temp[1, 2, ]
)
cor_dcc_prec <- rcor(dcc_prec)
cor_dcc_prec <- as.numeric(
  cor_dcc_prec[1, 2, ]
)


dcc_results <- data.frame(
  year = climate_rwi$year,
  Temperature = cor_dcc_temp,
  Precipitation = cor_dcc_prec
)

head(dcc_results)
tail(dcc_results)
summary(dcc_results[, c("Temperature","Precipitation")])
colSums(is.na(dcc_results))

## Figure DCC - RWI, TEMPERATURE E PRECIPITATION

library(ggplot2)

### Transformation into long format
dcc_long <- rbind(
  data.frame(
    year = dcc_results$year,
    correlation = dcc_results$Temperature,
    variable = "Temperature - RWI"
  ),
  data.frame(
    year = dcc_results$year,
    correlation = dcc_results$Precipitation,
    variable = "Precipitation - RWI"
  )
)



# GRAPH
p_dcc <- ggplot(
  dcc_long,
  aes(
    x = year,
    y = correlation,
    colour = variable
  )
) +
  geom_hline(
    yintercept = 0,
    colour = "black",
    linewidth = 0.5
  ) +
  geom_line(
    linewidth = 1.1
  ) +
  scale_colour_manual(
    values = c(
      "Temperature - RWI" = "#E63946",
      "Precipitation - RWI" = "#0077B6"
    ),
    name = NULL
  ) +
  scale_x_continuous(
    breaks = seq(1901, 2025, by = 10),
    expand = c(0, 0)
  ) +
  labs(
    title = "Dynamic Conditional Correlation (DCC)",
    subtitle = "RWI–Temperature and RWI–Precipitation",
    x = "Year",
    y = "Dynamic conditional correlation"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    
    axis.text.x = element_text(
      size = 8,
      colour = "black",
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(
      size = 9,
      colour = "black"
    ),
    axis.title.x = element_text(
      size = 10,
      face = "bold"
    ),
    axis.title.y = element_text(
      size = 10,
      face = "bold"
    ),
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      size = 10,
      hjust = 0.5
    ),
    legend.position = "bottom",
    legend.text = element_text(
      size = 10
    ),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.6
    ),
    plot.margin = margin(
      10,
      15,
      10,
      10
    )
  )

print(p_dcc)


# ALTERNATIVE DCC: LOG PRECIPITATION
### Logaritmo delle precipitazioni annuali
Precipitation_log <- log(
  climate_rwi$precipitation
)
### Dataset
data_dcc_prec_log <- cbind(
  RWI = climate_rwi$RWI,
  Precipitation = Precipitation_log
)
### DCC specification
spec_dcc_prec_log <- dccspec(
  uspec = multispec(
    replicate(
      2,
      spec_garch
    )
  ),
  dccOrder = c(1, 1),
  distribution = "mvnorm"
)

### Fit
dcc_prec_log <- dccfit(
  spec_dcc_prec_log,
  data = data_dcc_prec_log,
  fit.control = list(
    eval.se = TRUE
  )
)
show(dcc_prec_log)


## EXTRACT LOG-PRECIPITATION DCC
cor_dcc_prec_log <- rcor(
  dcc_prec_log
)
cor_dcc_prec_log <- as.numeric(
  cor_dcc_prec_log[1, 2, ]
)
summary(cor_dcc_prec_log)
range(cor_dcc_prec_log)
length(
  unique(
    round(
      cor_dcc_prec_log,
      6
    )
  )
)
