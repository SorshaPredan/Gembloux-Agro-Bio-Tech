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

# READING XCT DATA
#### make XCT.read function ####
XCT.read <- function(path,# A path to the folder containing the txt files
                     output = "ringwidth_density", # The output type, can be "ringwidth" (dplR format of ring width), "density" (dplR format of density parameter), "ringwidth_density" (long format of the sample, year, ring width, and density), or "density_profile" (long format of the sample, year, and density profile in that year)
                     densityType = "fraction", # The type of density to calculate, can be "fraction" or "fixed". "fraction" calculates the density in a variable width window that corresponds to two fraction numbers that go from 0 (start ring) to 1 (end ring), set in variable area. "fixed" calculates the density in a fixed width window, starting from the beginning or the end of the ring. set in variable area.
                     area = c(0.75, 1), # Fraction of the ring to calculate the density parameter. If densityType = "fraction" this is a vector of two numbers that go from 0 (start ring) to 1 (end ring). If densityType = "fixed" this is a vector with "start" or "end" as the first variable, and the width of the window in micrometers as the second variable.
                     fun = "mean", # The function to calculate the density in the selected area, can be "mean", "median", "min", "max", or "mean_top_x". "mean_top_x" calculates the mean of the x highest values in the selected area, the variable x should be set to a fraction between 0 and 1.
                     x = 0.2,  # Fraction of the highest values to calculate the mean. Only used if fun = "mean_top_x".
                     removeNarrowRings = FALSE, # Removes density parameters of rings that are too small, set in minRingWidth. Can be either TRUE or FALSE.
                     minRingWidth = 0.030, #  Minimum width of the ring in mm that should be used in density calculations, only if removeNarrowRings = TRUE
                     overruleResolution = FALSE, # Overrule the resolution of the XCT data txts. If TRUE, the resolution of the XCT data is set to the resolution parameter. If FALSE, the resolution is set to the value in the ringwidth.txt file.
                     resolution = 1 # The resolution of the data in µm/pixel. Only used if overruleResolution = TRUE.
                     ){

   
  # Load all ringwidth.txt files
  files <- list.files(path, pattern = "_ringwidth.txt", full.names = TRUE)
  # Function to read each file, add the file name, and combine all into one data frame
  rings <- lapply(files, function(file) {
    # Extract the file name before "_density_corr"
    file_name <- sub("_ringwidth.txt$", "", basename(file))
    # Read the file, treat "NaN" as NA, and specify numeric columns
    data <- read_delim(file, delim = ", ", 
                       col_names = c("width", "Year", "pixelsize", "Felldate", "MissingringsBefore", "BrokenRingType"),
                       na = "NaN", show_col_types = FALSE)  
    data <- data %>% mutate(Sample = file_name, row_number = row_number()) # Add the file name as a new column
    data <- data %>% filter(BrokenRingType != 2) # remove type 2 broken rings rows
    data <- data %>% group_by(Year) %>% mutate(width = max(width)) %>% ungroup() # set ringwidth to the max ringwidth: relevant for ring a type 1 broken ring present
    data <- data %>% filter(BrokenRingType != 1) # remove type 1 broken rings rows
    if (overruleResolution) {data$pixelsize <- resolution} # set the resolution to the value in the resolution parameter
    data$RW <- data$width * data$pixelsize /1000 # calculate ring width in mm
    return(data)
  }) %>% bind_rows() # Combine all data frames into one
  
  if (output == "ringwidth") {
    # Return the ring width data frame in dplR format
    rings <- rings[, c("Year", "Sample", "RW")]
    rings <- rings %>% na.omit()
    rings <- rings %>% group_by(Sample, Year) %>% summarise(RW = max(RW), .groups = "drop") # merge duplicate rows (relevant for type 1 broken rings)
    rings <- pivot_wider(rings, names_from = "Sample", values_from = "RW")
    rings <- rings %>% arrange(Year)
    rings <- rings %>% complete(Year = seq(min(rings$Year), max(rings$Year), 1)) # insert NA column in years that are skipped so there is a column for each year
    extent <- as.vector(rings$Year) # year extent data
    rings$Year <- NULL # delete year column
    rings <- as.data.frame(rings) # change to dataframe
    row.names(rings) <- extent # change row names to years
    return(rings)
  }
  
  # Load all _density_corr.txt files
  files <- list.files(path, pattern = "_density_corr.txt", full.names = TRUE)
  # Function to read each file, add the file name, and combine all into one data frame
  Density_corr <- lapply(files, function(file) {
    # Extract the file name before "_density_corr"
    file_name <- sub("_density_corr.txt$", "", basename(file)) 
    # Read the file, treat "NaN" as NA
    data <- read_delim(file, delim = "\n", col_names = "Density", na = "NaN",  show_col_types = FALSE)
    # Add the file name as a new column
    data <- data %>% mutate(Sample = file_name, row_number = row_number())
    return(data)
  }) %>% bind_rows() # Combine all data frames into one
  
  # Load all _zpos_corr files
  files <- list.files(path, pattern = "_zpos_corr.txt", full.names = TRUE)
  # Function to read each file, add the file name, and combine all into one data frame
  zpos_corr <- lapply(files, function(file) {
    # Extract the file name before "_density_corr"
    file_name <- sub("_zpos_corr.txt$", "", basename(file))
    # Read the file, treat "NaN" as NA
    data <- read_delim(file, delim = "\n", col_names = "xpos", na = "NaN",  show_col_types = FALSE)
    # Add the file name as a new column
    data <- data %>% mutate(Sample = file_name, row_number = row_number()-1)
    return(data)
  }) %>% bind_rows() # Combine all data frames into one

  # Merge rings and zpos_corr by Sample and row_number to add start and end columns
  rings <- rings %>%
    left_join(zpos_corr %>% 
                mutate(row_number = row_number + 1) %>% # Shift row_number for "start"
                rename(start = xpos), 
              by = c("Sample", "row_number")) %>%
    left_join(zpos_corr %>% rename(end = xpos), by = c("Sample", "row_number"))
  rings$start <- rings$start + 1 # first pixel is the pixel after the end of the last one, otherwise this pixel is accounted twice
  rings <- rings %>% na.omit()   # remove NAs
  ringsbackup <- rings
  
  if (removeNarrowRings) {
    rings <- rings %>% filter(RW >= minRingWidth)  # remove rings with width smaller than minRingWidth
  }
  
  # put year and Sample in density_corr
  density_map <- rings %>% # Create a helper data frame to link `Density_corr` with `rings`, This creates a mapping of Sample, Year, and the range of row_numbers for each ring
    select(Sample, Year, start, end) %>%
    distinct() %>%
    rowwise() %>%
    mutate(row_number = list(seq(start, end))) %>%  # Generate a sequence for each range
    unnest(cols = c(row_number))  # Expand each sequence to individual rows
  Density_corr <- Density_corr %>%   # Join this mapping to `Density_corr` by Sample and row_number
    left_join(density_map, by = c("Sample", "row_number")) %>%
    arrange(Sample, row_number) %>%
    group_by(Sample, Year) %>%
    mutate(row_number = row_number()) %>% # Recalculate row_number continuously within each ring
    ungroup()
  Density_corr <- Density_corr[!is.na(Density_corr$Year),]   # remove NAs (gaps in density profile)
  
  if (output == "density_profile") {
    # Return the density profile data frame in long format
    Density_corr <- Density_corr[, c("Sample", "Year", "row_number", "Density")] %>% 
      rename(Pixel_nr_along_ring = row_number)   %>%  
      arrange(Sample, Year, Pixel_nr_along_ring) %>%
      group_by(Sample) %>% mutate(row_number_along_sample = row_number()) %>% ungroup()    # extra column that gives the pixel number along one whole sample
    return(Density_corr)
  }
  
  # Function to calculate the mean of the top x values in a vector
  mean_top_x <- function(vec, x) {
    # Ensure x is between 0 and 100
    if (x < 0 || x > 1) {
      stop("x should be a fraction between 0 and 1")}
    # Remove any NA values
    vec <- na.omit(vec)
    # Calculate the number of top values to consider
    n_top <- ceiling(length(vec) * x)
    # Sort the vector in descending order and take the top n values
    top_values <- sort(vec, decreasing = TRUE)[1:n_top]
    # Calculate and return the mean of the top values
    return(ifelse( !all(is.na(top_values)), mean(top_values, na.rm=TRUE), NA))
  }
  # Custom function for calculating the density based on `fun` parameter
  calculate_density <- function(density_values, fun, x) {
    
    if (fun == "mean") {
      return(ifelse( !all(is.na(density_values)), mean(density_values, na.rm=TRUE), NA))
    } else if (fun == "median") {
      return(ifelse( !all(is.na(density_values)), median(density_values, na.rm=TRUE), NA))
    } else if (fun == "min") {
      return(ifelse( !all(is.na(density_values)), min(density_values, na.rm=TRUE), NA))
    } else if (fun == "max") {
      return(ifelse( !all(is.na(density_values)), max(density_values, na.rm=TRUE), NA))
    } else if (fun == "mean_top_x") {
      # Use mean_top_x function for the top x% of values
      return(mean_top_x(density_values, x))
    } else {
      stop("Invalid function specified in `fun` argument.")
    }
  }
  
  # calculate density of fraction
  if (densityType == "fraction") {
    Density_corr <- Density_corr %>%
      group_by(Sample, Year) %>%
      summarise(
        Density = calculate_density(Density[seq_along(Density) > ((area[1]) * length(Density)) & seq_along(Density) <= ((area[2]) * length(Density))], fun = fun, x = x),
        .groups = "drop"  # Ungroups the output completely
      )
  }
  
  # calculate density in fixed area
  if (densityType == "fixed") {
    start_or_end <- area[1]
    lengthMicron <- as.numeric(area[2])
    if (start_or_end == "start") {
      Density_corr <- Density_corr %>%
        group_by(Sample, Year) %>%
        summarise(
          Density = calculate_density(Density[seq_along(Density) <= round(lengthMicron / (mean(rings$pixelsize)))], fun = fun, x = x),
          .groups = "drop"  # Ungroups the output completely
        )
    }
    if (start_or_end == "end") {
      Density_corr <- Density_corr %>%
        group_by(Sample, Year) %>%
        summarise(
          Density = calculate_density(Density[seq_along(Density) > (length(Density) - round(lengthMicron / (mean(rings$pixelsize))))], fun = fun, x = x),
          .groups = "drop"  # Ungroups the output completely
        )
    }
  }
  if (output == "density") {
    # Return the density data frame in dplR format
    Density_corr <- Density_corr[, c("Year", "Sample", "Density")]
    Density_corr <- Density_corr %>% na.omit()
    Density_corr <- pivot_wider(Density_corr, names_from = "Sample", values_from = "Density")
    Density_corr <- Density_corr %>% arrange(Year)
    Density_corr <- Density_corr %>% complete(Year = seq(min(Density_corr$Year), max(Density_corr$Year), 1)) # insert NA column in years that are skipped so there is a column for each year
    extent <- as.vector(Density_corr$Year) # year extent data
    Density_corr$Year <- NULL # delete year column
    Density_corr <- as.data.frame(Density_corr) # change to dataframe
    row.names(Density_corr) <- extent # change row names to years
    return(Density_corr)
  }
  if (output == "ringwidth_density") {
    # Return the ring width and density data frame in long format
    rings <- ringsbackup  %>%  group_by(Sample, Year) %>% summarise(RW = max(RW), .groups = "drop")  # merge duplicate rows (relevant for type 1 broken rings)
    Data <- Density_corr %>% left_join(rings %>% select(Sample, Year, RW), by = c("Sample", "Year"))
    return(Data)
  }  
}

### Which months most strongly influence EWD, LWD, and MXD? ###
## PC Gembloux
XCT_folder <- "C:/Users/user/Desktop/TIROCINIO FINALE/DENSITA/XCT"
## PC Sorsha
XCT_folder <- "C:/Users/user/Desktop/TIROCINIO FINALE/DENSITA/XCT/CORES"
# EARLYWOOD DENSITY
EWD <- XCT.read(
  path = XCT_folder,
  output = "density",
  densityType = "fraction",
  area = c(0,0.25),
  fun = "mean"
)
head(EWD)
dim(EWD)
class(EWD)
class(EWD) <- c("rwl","data.frame")
EWDChron <- chron(EWD)
## EWD + Prec/Temp
EWD_Prec <- dcc(EWDChron, PrecSites, selection = -6:9, method = "correlation",
                   timespan = c(1930,1990), var_names = "Precipitation", boot = "std")
class(EWD_Prec)
plot(EWD_Prec)
EWD_Temp <- dcc(EWDChron, TempSites, selection = -6:9, method = "correlation",
                   timespan = c(1930,1990), var_names = "Temperature", boot = "std")
class(EWD_Temp)
plot(EWD_Temp)
## EWD correlation plot
prec <- EWD_Prec$coef$coef
temp <- EWD_Temp$coef$coef
mesi <- EWD_Prec$coef$month
### matrice correlazioni
corr_matrix <- rbind(
  prec,
  temp
)
### grafico
bar_position <- barplot(
  corr_matrix,
  beside = TRUE,
  names.arg = mesi,
  col = c("steelblue", "red"),
  ylim = c(-0.75,0.75),
  ylab = "Pearson r",
  xlab = "Month",
  las = 2,
  main = "EWD - Climate correlations"
)
abline(h = 0, lwd = 2)
legend(
  "topright",
  legend = c("Precipitation", "Temperature"),
  fill = c("steelblue", "red"),
  bty = "n"
)
sig_prec <- EWD_Prec$coef$significant
sig_temp <- EWD_Temp$coef$significant
pos_prec <- bar_position[1,]
pos_temp <- bar_position[2,]
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

# LATEWOOD DENSITY
LWD <- XCT.read(
  path = XCT_folder,
  output = "density",
  densityType = "fraction",
  area = c(0.75,1),
  fun = "mean"
)
head(LWD)
dim(LWD)
class(LWD)
class(LWD) <- c("rwl","data.frame")
LWDChron <- chron(LWD)
## LWD + Prec/Temp
LWD_Prec <- dcc(LWDChron, PrecSites, selection = -6:9, method = "correlation",
                   timespan = c(1930,1990), var_names = "Precipitation", boot = "std")
class(LWD_Prec)
plot(LWD_Prec)
LWD_Temp <- dcc(LWDChron, TempSites, selection = -6:9, method = "correlation",
                   timespan = c(1930,1990), var_names = "Temperature", boot = "std")
class(LWD_Temp)
plot(LWD_Temp)
## LWD correlation plot
prec <- LWD_Prec$coef$coef
temp <- LWD_Temp$coef$coef
mesi <- LWD_Prec$coef$month
### matrice correlazioni
corr_matrix <- rbind(
  prec,
  temp
)
### grafico
bar_position <- barplot(
  corr_matrix,
  beside = TRUE,
  names.arg = mesi,
  col = c("steelblue", "red"),
  ylim = c(-0.75,0.75),
  ylab = "Pearson r",
  xlab = "Month",
  las = 2,
  main = "LWD - Climate correlations"
)
abline(h = 0, lwd = 2)
legend(
  "topright",
  legend = c("Precipitation", "Temperature"),
  fill = c("steelblue", "red"),
  bty = "n"
)
sig_prec <- LWD_Prec$coef$significant
sig_temp <- LWD_Temp$coef$significant
pos_prec <- bar_position[1,]
pos_temp <- bar_position[2,]
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

# MAXIMUM DENSITY
MXD <- XCT.read(
  path = XCT_folder,
  output = "density",
  densityType = "fraction",
  area = c(0,1),
  fun = "max"
)
head(MXD)
dim(MXD)
class(MXD)
class(MXD) <- c("rwl","data.frame")
MXDChron <- chron(MXD)
## MXD + Prec/Temp
MXD_Prec <- dcc(MXDChron, PrecSites, selection = -6:9, method = "correlation",
                   timespan = c(1930,1990), var_names = "Precipitation", boot = "std")
class(MXD_Prec)
plot(MXD_Prec)
MXD_Temp <- dcc(MXDChron, TempSites, selection = -6:9, method = "correlation",
                   timespan = c(1901,2025), var_names = "Temperature", boot = "std", dynamic="movying")
class(MXD_Temp)
plot(MXD_Temp)

## MXD correlation plot
prec <- MXD_Prec$coef$coef
temp <- MXD_Temp$coef$coef
mesi <- MXD_Prec$coef$month
### matrice correlazioni
corr_matrix <- rbind(
  prec,
  temp
)
### grafico
bar_position <- barplot(
  corr_matrix,
  beside = TRUE,
  names.arg = mesi,
  col = c("steelblue", "red"),
  ylim = c(-0.75,0.75),
  ylab = "Pearson r",
  xlab = "Month",
  las = 2,
  main = "MXD - Climate correlations"
)
abline(h = 0, lwd = 2)
legend(
  "topright",
  legend = c("Precipitation", "Temperature"),
  fill = c("steelblue", "red"),
  bty = "n"
)
sig_prec <- MXD_Prec$coef$significant
sig_temp <- MXD_Temp$coef$significant
pos_prec <- bar_position[1,]
pos_temp <- bar_position[2,]
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

# PRECIPITATION AND TEMPERATURE CHRONOLOGIES
## Conversion of chronologies into density data frames
EWD_df <- data.frame(
  year = as.numeric(rownames(EWDChron)),
  density = EWDChron$std,
  parameter = "EWD"
)
LWD_df <- data.frame(
  year = as.numeric(rownames(LWDChron)),
  density = LWDChron$std,
  parameter = "LWD"
)
MXD_df <- data.frame(
  year = as.numeric(rownames(MXDChron)),
  density = MXDChron$std,
  parameter = "MXD"
)
density_all <- bind_rows(EWD_df, LWD_df, MXD_df)
head(density_all)
table(density_all$parameter)

## Wood density chronologies over time
library(ggplot2)
p_density <- ggplot(
  density_all,
  aes(
    x = year,
    y = density,
    colour = parameter
  )
) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(
    values = c(
      EWD = "darkgreen",
      LWD = "orange",
      MXD = "red"
    )
  ) +
  labs(
    x = "Year",
    y = "Density",
    colour = "Chronology",
    title = "Wood density chronologies"
  ) +
  theme_classic()
p_density

# Precipitation and Temperature plot 
## Annual precipitation
Prec_year <- data.frame(
  year = PrecSites$year,
  precipitation = rowSums(PrecSites[,2:13], na.rm = TRUE)
)
head(Prec_year)
dim(Prec_year)
## Annual temperature
Temp_year <- data.frame(
  year = TempSites$year,
  temperature = rowMeans(TempSites[,2:13], na.rm = TRUE)
)
head(Temp_year)
dim(Temp_year)
### Dual-axis plot
library(patchwork)
Climate_year <- merge(
  Prec_year,
  Temp_year,
  by="year"
)
head(Climate_year)
scale_factor <- max(Climate_year$precipitation) /
                max(Climate_year$temperature)
Climate_year$temp_scaled <- Climate_year$temperature * scale_factor
p_climate <- ggplot(Climate_year, aes(x=year)) +
  geom_line(
    aes(y=precipitation),
    colour="steelblue",
    linewidth=0.8
  ) +
  geom_line(
    aes(y=temp_scaled),
    colour="red",
    linewidth=0.8
  ) +
  scale_y_continuous(
    name="Precipitation",
    sec.axis = sec_axis(
      ~./scale_factor,
      name="Temperature"
    )
  ) +
  labs(
    x="Year"
  ) +
  theme_classic()
p_climate

p_density / p_climate

### Has the relationship between temperature, precipitation, and EWD/LWD/MXD remained stable over time, or has it changed? ###
## Summer temperature (JJA)
Temp_JJA <- data.frame(
  year = TempSites$year,
  temperature = rowMeans(
    TempSites[,c("Jun","Jul","Aug")],
    na.rm = TRUE
  )
)
head(Temp_JJA)
Temp_JJA$z_temp <- as.numeric(scale(Temp_JJA$temperature))

### EWD Z-SCORE
EWD_z <- data.frame(
  year = as.numeric(rownames(EWDChron)),
  z = as.numeric(scale(EWDChron$std)),
  parameter = "EWD"
)
head(EWD_z)
### EWD + TEMPERATURA
EWD_Temp_JJA <- merge(
  EWD_z,
  Temp_JJA,
  by = "year"
)
head(EWD_Temp_JJA)

### LWD Z-SCORE
LWD_z <- data.frame(
  year = as.numeric(rownames(LWDChron)),
  z = as.numeric(scale(LWDChron$std)),
  parameter = "LWD"
)
head(LWD_z)
### LWD + TEMPERATURA
LWD_Temp_JJA <- merge(
  LWD_z,
  Temp_JJA,
  by = "year"
)

### MXD Z-SCORE
MXD_z <- data.frame(
  year = as.numeric(rownames(MXDChron)),
  z = as.numeric(scale(MXDChron$std)),
  parameter = "MXD"
)
head(MXD_z)
### MXD + TEMPERATURA
MXD_Temp_JJA <- merge(
  MXD_z,
  Temp_JJA,
  by = "year"
)
head(MXD_Temp_JJA)

# Z-score correlation analysis
## EWD - temperature JJA
EWD_cor <- cor.test(
  EWD_Temp_JJA$z,
  EWD_Temp_JJA$z_temp,
  method = "pearson"
)
## LWD - temperature JJA
LWD_cor <- cor.test(
  LWD_Temp_JJA$z,
  LWD_Temp_JJA$z_temp,
  method = "pearson"
)
## MXD - temperature JJA
MXD_cor <- cor.test(
  MXD_Temp_JJA$z,
  MXD_Temp_JJA$z_temp,
  method = "pearson"
)
### Value extraction
EWD_r <- EWD_cor$estimate
EWD_p <- EWD_cor$p.value
LWD_r <- LWD_cor$estimate
LWD_p <- LWD_cor$p.value
MXD_r <- MXD_cor$estimate
MXD_p <- MXD_cor$p.value
EWD_r
EWD_p
LWD_r
LWD_p
MXD_r
MXD_p

# Plot EWD - Temperature JJA ###
library(ggplot2)
ggplot(EWD_Temp_JJA, aes(x = year)) +
  ### EWD Z-score variability band
  geom_ribbon(
    aes(
      ymin = z - sd(z),
      ymax = z + sd(z)
    ),
    fill = "darkgreen",
    alpha = 0.15
  ) +
  ### line EWD
  geom_line(
    aes(y = z),
    colour = "darkgreen",
    linewidth = 1
  ) +
  ### line temperature
  geom_line(
    aes(y = z_temp),
    colour = "red",
    linewidth = 1
  ) +
  ### line at zero
  geom_hline(
    yintercept = 0,
    colour = "grey50",
    linetype = "dashed"
  ) +
  scale_y_continuous(
    limits = c(-2.75,2.75)
  ) +
  labs(
    title = paste0(
      "EWD - Summer temperature relationship\n",
      "Pearson r = ",
      round(EWD_r,2),
      "   p = ",
      round(EWD_p,3)
    ),
    x = "Year",
    y = "Z-score",
    caption = "Green: EWD chronology | Red: summer temperature (JJA)"
  ) +
  theme_classic()

### Plot LWD - Temperature JJA ###
ggplot(LWD_Temp_JJA, aes(x = year)) +
  geom_ribbon(
    aes(
      ymin = z - sd(z),
      ymax = z + sd(z)
    ),
    fill = "orange",
    alpha = 0.15
  ) +
  geom_line(
    aes(y = z),
    colour = "orange",
    linewidth = 1
  ) +
  geom_line(
    aes(y = z_temp),
    colour = "red",
    linewidth = 1
  ) +
  geom_hline(
    yintercept = 0,
    colour = "grey50",
    linetype = "dashed"
  ) +
  scale_y_continuous(
    limits = c(-2.75,2.75)
  ) +
  labs(
    title = paste0(
      "LWD - Summer temperature relationship\n",
      "Pearson r = ",
      round(LWD_r,2),
      "   p = ",
      round(LWD_p,3)
    ),
    x = "Year",
    y = "Z-score",
    caption = "Orange: LWD chronology | Red: summer temperature (JJA)"
  ) +
  theme_classic()
