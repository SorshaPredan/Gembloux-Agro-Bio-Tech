# TRUE HEATMAP - 25 YEAR MOVING CORRELATION
# SEASONAL DATA WITH SPECIFIED MONTHS

library(ggplot2)

ordine <- c(
  "temp.mean.curr.winter",
  "temp.mean.curr.spring",
  "temp.mean.curr.summer",
  "prec.curr.winter",
  "prec.curr.spring",
  "prec.curr.summer"
)
grafico$variabile <- factor(
  grafico$variabile,
  levels = rev(ordine)
)

## labels with months 
etichette <- c(
  "temp.mean.curr.winter" =
    "Temp mean curr. Winter (Dec-Jan-Feb)",
  "temp.mean.curr.spring" =
    "Temp mean curr. Spring (Mar-Apr-May)",
  "temp.mean.curr.summer" =
    "Temp mean curr. Summer (Jun-Jul-Aug)",
  "prec.curr.winter" =
    "Prec curr. Winter (Dec-Jan-Feb)",
  "prec.curr.spring" =
    "Prec curr. Spring (Mar-Apr-May)",
  "prec.curr.summer" =
    "Prec curr. Summer (Jun-Jul-Aug)"
)

## correlation classes 
grafico$classe <- cut(
  grafico$coefficiente,
  breaks = c(
    -1,
    0,
    0.50,
    0.70,
    0.80,
    0.85,
    0.90,
    0.95,
    1
  ),
  labels = c(
    "< 0",
    "0-0.50",
    "0.50-0.70",
    "0.70-0.80",
    "0.80-0.85",
    "0.85-0.90",
    "0.90-0.95",
    "0.95-1.00"
  ),
  include.lowest = TRUE
)

## heatmap
p <- ggplot(
  grafico,
  aes(
    x = anno,
    y = variabile,
    fill = classe
  )
) 

## Heatmap tiles
  geom_tile(
    width = 1,
    height = 0.95,
    colour = "white",
    linewidth = 0.25
  ) +

## colors 
  scale_fill_manual(
    values = c(
      "< 0"       = "#2166AC",
      "0-0.50"    = "#67A9CF",
      "0.50-0.70" = "#D1E5F0",
      "0.70-0.80" = "#F7F7F7",
      "0.80-0.85" = "#FDDDA0",
      "0.85-0.90" = "#F4A582",
      "0.90-0.95" = "#D6604D",
      "0.95-1.00" = "#B2182B"
    ),
    na.value = "#BDBDBD",
    drop = FALSE,
    name = "Coeff."
  ) +

## x axis 5. X AXIS
  scale_x_continuous(
    breaks = seq(
      1901,
      1990,
      by = 5
    ),
    expand = c(0, 0)
  ) +

## y axis 6. Y AXIS
  scale_y_discrete(
    labels = etichette
  ) +

## titles
  labs(
    title = "25-Year Moving Correlation",
    x = "Time",
    y = NULL
  ) +

## theme 
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(
      size = 9,
      colour = "black"
    ),
    axis.text.x = element_text(
      size = 8,
      colour = "black",
      angle = 45,
      hjust = 1
    ),
    axis.title.x = element_text(
      size = 10,
      face = "bold"
    ),
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5
    ),
    legend.title = element_text(
      size = 10,
      face = "bold"
    ),
    legend.text = element_text(
      size = 8
    ),
    panel.border = element_blank(),
    plot.margin = margin(
      10,
      15,
      10,
      10
    )
  )

print(p)
