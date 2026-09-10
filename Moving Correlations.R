

############################################################
# 4. EXPORT RING-WIDTH DATA
############################################################

write.xlsx(
  TRW,
  "TRW.xlsx",
  rowNames = TRUE
)


############################################################
# 5. BASAL AREA INCREMENT
############################################################

BeechBAI <- bai.in(TRW)


############################################################
# 6. RAW TREE-RING SERIES
############################################################

matplot(
  as.numeric(rownames(TRW)),
  TRW,
  type = "l",
  xlab = "Year",
  ylab = "Tree-ring width",
  main = "Raw tree-ring series"
)


############################################################
# 7. STANDARDIZATION
############################################################

TRWdetrend <- detrend(
  TRW,
  method = "Spline",
  nyrs = 30
)


############################################################
# 8. CHRONOLOGY
############################################################

BeechChron <- chron(
  TRWdetrend,
  prefix = "AVG",
  biweight = TRUE,
  prewhiten = FALSE
)

# Plot chronology
plot.crn(
  BeechChron,
  main = "Beech chronology"
)

# Check chronology time span
range(time(BeechChron))


############################################################
# 9. PREPARATION OF PRECIPITATION DATA
############################################################

# Set column names
colnames(PrecSites) <- c(
  "year",
  "Jan", "Feb", "Mar", "Apr", "May", "Jun",
  "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
)

# Convert monthly precipitation values to numeric
PrecSites[, 2:13] <- lapply(
  PrecSites[, 2:13],
  function(x) as.numeric(gsub("\\.", "", x))
)

# Apply original dataset scaling
PrecSites[, 2:13] <- PrecSites[, 2:13] / 1000000

# Check for duplicated years
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

# Remove duplicated years
PrecSites <- PrecSites[
  !duplicated(PrecSites$year),
]


############################################################
# 10. PREPARATION OF TEMPERATURE DATA
############################################################

# Set column names
colnames(TempSites) <- c(
  "year",
  "Jan", "Feb", "Mar", "Apr", "May", "Jun",
  "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
)

# Convert monthly temperature values to numeric
TempSites[, 2:13] <- lapply(
  TempSites[, 2:13],
  function(x) as.numeric(gsub("\\.", "", x))
)

# Apply original dataset scaling
TempSites[, 2:13] <- TempSites[, 2:13] / 1000000

# Check for duplicated years
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

# Remove duplicated years
TempSites <- TempSites[
  !duplicated(TempSites$year),
]


############################################################
# 11. DENDROCLIMATIC ANALYSIS
#
# DCC = Dendro Climatic Analysis
#
# dynamic = "moving" performs the moving-window analysis.
############################################################

# Climate selection:
# -6 = June of previous year
# -5 = July of previous year
# -4 = August of previous year
# -3 = September of previous year
# -2 = October of previous year
# -1 = November of previous year
#  1 = January of current year
#  2 = February of current year
# ...
#  9 = September of current year

climate_selection <- -6:9


############################################################
# 12. MOVING DCC - PRECIPITATION
############################################################

PrecCorr <- dcc(
  BeechChron,
  PrecSites,
  selection = climate_selection,
  method = "correlation",
  dynamic = "moving",
  timespan = c(1930, 1990),
  var_names = "precipitation",
  boot = "std"
)

# Plot provided by treeclim
plot(
  PrecCorr,
  main = "Moving climate-growth correlation - Precipitation"
)


############################################################
# 13. MOVING DCC - TEMPERATURE
############################################################

TempCorr <- dcc(
  BeechChron,
  TempSites,
  selection = climate_selection,
  method = "correlation",
  dynamic = "moving",
  timespan = c(1930, 1990),
  var_names = "temperature",
  boot = "std"
)

# Plot provided by treeclim
plot(
  TempCorr,
  main = "Moving climate-growth correlation - Temperature"
)


############################################################
# 14. EXTRACT MOVING CORRELATION COEFFICIENTS
############################################################

# Extract correlation coefficients
prec_coef <- PrecCorr$coef$coef
temp_coef <- TempCorr$coef$coef

# Extract significance matrices
prec_sig <- PrecCorr$coef$significant
temp_sig <- TempCorr$coef$significant




############################################################
# 15. PREPARING DATA FOR THE HEATMAP
############################################################

# Extract correlation matrices from treeclim
prec_matrix <- PrecCorr$coef$coef
temp_matrix <- TempCorr$coef$coef

# Extract significance matrices
prec_sig_matrix <- PrecCorr$coef$significant
temp_sig_matrix <- TempCorr$coef$significant


############################################################
# 16. CONVERT CORRELATION MATRICES TO DATA FRAMES
############################################################

# Precipitation
prec_heat <- as.data.frame(prec_matrix)

prec_heat$climate_month <- rownames(prec_matrix)

prec_heat <- prec_heat %>%
  pivot_longer(
    cols = -climate_month,
    names_to = "moving_window",
    values_to = "correlation"
  )


# Temperature
temp_heat <- as.data.frame(temp_matrix)

temp_heat$climate_month <- rownames(temp_matrix)

temp_heat <- temp_heat %>%
  pivot_longer(
    cols = -climate_month,
    names_to = "moving_window",
    values_to = "correlation"
  )


############################################################
# 17. CONVERT SIGNIFICANCE MATRICES TO DATA FRAMES
############################################################

# Precipitation
prec_sig <- as.data.frame(prec_sig_matrix)

prec_sig$climate_month <- rownames(prec_sig_matrix)

prec_sig <- prec_sig %>%
  pivot_longer(
    cols = -climate_month,
    names_to = "moving_window",
    values_to = "significant"
  )


# Temperature
temp_sig <- as.data.frame(temp_sig_matrix)

temp_sig$climate_month <- rownames(temp_sig_matrix)

temp_sig <- temp_sig %>%
  pivot_longer(
    cols = -climate_month,
    names_to = "moving_window",
    values_to = "significant"
  )


############################################################
# 18. COMBINE CORRELATION AND SIGNIFICANCE DATA
############################################################

# Precipitation
prec_heat <- prec_heat %>%
  left_join(
    prec_sig,
    by = c(
      "climate_month",
      "moving_window"
    )
  )


# Temperature
temp_heat <- temp_heat %>%
  left_join(
    temp_sig,
    by = c(
      "climate_month",
      "moving_window"
    )
  )


############################################################
# 19. DEFINE CLIMATE MONTH ORDER
############################################################

# Precipitation
month_order <- c(
  "precipitation.prev.jun",
  "precipitation.prev.jul",
  "precipitation.prev.aug",
  "precipitation.prev.sep",
  "precipitation.prev.oct",
  "precipitation.prev.nov",
  "precipitation.prev.dec",
  "precipitation.jan",
  "precipitation.feb",
  "precipitation.mar",
  "precipitation.apr",
  "precipitation.may",
  "precipitation.jun",
  "precipitation.jul",
  "precipitation.aug",
  "precipitation.sep"
)


# Temperature
temp_month_order <- c(
  "temperature.prev.jun",
  "temperature.prev.jul",
  "temperature.prev.aug",
  "temperature.prev.sep",
  "temperature.prev.oct",
  "temperature.prev.nov",
  "temperature.prev.dec",
  "temperature.jan",
  "temperature.feb",
  "temperature.mar",
  "temperature.apr",
  "temperature.may",
  "temperature.jun",
  "temperature.jul",
  "temperature.aug",
  "temperature.sep"
)






############################################################
# 20. HEATMAP - PRECIPITATION
############################################################

# Set the order of climate months
prec_heat$climate_month <- factor(
  prec_heat$climate_month,
  levels = month_order
)

# Create precipitation heatmap
heatmap_prec <- ggplot(
  prec_heat,
  aes(
    x = moving_window,
    y = climate_month,
    fill = correlation
  )
) +

  # Heatmap tiles
  geom_tile(
    color = "#4D4D4D",
    linewidth = 0.35
  ) +

  # Mark significant correlations
  geom_text(
    aes(
      label = ifelse(significant, "*", "")
    ),
    color = "black",
    size = 3
  ) +

  # Correlation color scale
  scale_fill_gradient2(
    low = "#0047AB",
    mid = "#FFFFFF",
    high = "#D00000",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Pearson r"
  ) +

  # Climate month labels
  scale_y_discrete(
    labels = c(
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
  ) +

  # Labels
  labs(
    title = "Moving climate-growth correlations",
    subtitle = "Precipitation",
    x = "Moving 25-year window",
    y = "Climate month"
  ) +

  # Theme
  theme_minimal(base_size = 12) +

  theme(
    panel.grid = element_blank(),

    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      size = 8
    ),

    axis.text.y = element_text(
      size = 9
    ),

    plot.title = element_text(
      face = "bold"
    ),

    plot.subtitle = element_text(
      face = "italic"
    ),

    legend.position = "right"
  )

# Display precipitation heatmap
print(heatmap_prec)


############################################################
# 21. HEATMAP - TEMPERATURE
############################################################

# Set the order of climate months
temp_heat$climate_month <- factor(
  temp_heat$climate_month,
  levels = temp_month_order
)

# Create temperature heatmap
heatmap_temp <- ggplot(
  temp_heat,
  aes(
    x = moving_window,
    y = climate_month,
    fill = correlation
  )
) +

  # Heatmap tiles
  geom_tile(
    color = "#4D4D4D",
    linewidth = 0.35
  ) +

  # Mark significant correlations
  geom_text(
    aes(
      label = ifelse(significant, "*", "")
    ),
    color = "black",
    size = 3
  ) +

  # Correlation color scale
  scale_fill_gradient2(
    low = "#0047AB",
    mid = "#FFFFFF",
    high = "#D00000",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Pearson r"
  ) +

  # Climate month labels
  scale_y_discrete(
    labels = c(
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
  ) +

  # Labels
  labs(
    title = "Moving climate-growth correlations",
    subtitle = "Temperature",
    x = "Moving 25-year window",
    y = "Climate month"
  ) +

  # Theme
  theme_minimal(base_size = 12) +

  theme(
    panel.grid = element_blank(),

    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      size = 8
    ),

    axis.text.y = element_text(
      size = 9
    ),

    plot.title = element_text(
      face = "bold"
    ),

    plot.subtitle = element_text(
      face = "italic"
    ),

    legend.position = "right"
  )

# Display temperature heatmap
print(heatmap_temp)
