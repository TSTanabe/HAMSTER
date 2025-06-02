#!/usr/bin/env Rscript

# plot_neighborhood_confusion_A3_aligned_fixedMargin_filtered.R
#
# Dieses Skript splittet einen langen Neighborhood‐Confusion‐Plot in mehrere
# DIN A3‐Quer‐Seiten auf, filtert jedoch alle Neighborhoods heraus,
# die weder TP noch FN enthalten (TP == 0 und FN == 0).
# Änderungen gegenüber dem Originalskript:
#   • Filter: Nur Neighborhoods mit TP > 0 oder FN > 0 bleiben erhalten.
#   • Keine Balken‐Umrandung (kein color="black").
#   • Bar‐Länge wird auf 5000 begrenzt → x‐Achse fix [0,5000].
#   • Werte > 5000 werden fürs Label‐Placement (PlotY) auf 5000 gecapped.
#   • Abstand zwischen Balken = 0,15 cm (früher 0,10 cm), Balkenhöhe nach wie vor 0,25 cm.
#   • Der Rest bleibt: 12 pt‐Schrift, fixe linke Margin = 40 mm, A3‐Quer, Pagination.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(grid)
})

# -----------------------------------------------------------
# 1) Kommandozeile: Pfad zur TSV-Datei
# -----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Please provide the path to the TSV file as the first argument.\n  Example: Rscript plot_neighborhood_confusion.R /path/to/file.tsv")
}
tsv_path <- args[1]

# Default: enable filtering unless args[2] exists and is 'FALSE'
if (length(args) >= 2) {
  # as.logical("TRUE") → TRUE, as.logical("FALSE") → FALSE
  filter_flag <- as.logical(args[2])
  if (is.na(filter_flag)) {
    stop("The second argument must be 'TRUE' or 'FALSE', but '", args[2], "' was provided.")
  }
} else {
  filter_flag <- TRUE
}

# -----------------------------------------------------------
# 2) Output-Pfad berechnen und prüfen, ob bereits vorhanden
# -----------------------------------------------------------
# Bestimme das Suffix: "_reference_csb.pdf" wenn filter_flag == TRUE, sonst "_all_csb.pdf"
suffix <- if (filter_flag) "_reference_csb.pdf" else "_all_csb.pdf"
# Ersetze die .tsv-Endung durch das jeweilige Suffix
output_file <- sub("\\.tsv$", suffix, tsv_path)

if (file.exists(output_file)) {
  message("[INFO] Output file '", output_file, "' already exists – skipping generation.")
  quit(status = 0)
}

# -----------------------------------------------------------
# 3) Lese TSV: Spalten = Neighborhood, TP, FP, FN, TN, Total
# -----------------------------------------------------------
df <- read_tsv(
  tsv_path,
  col_types = cols(
    Neighborhood = col_character(),
    TP           = col_double(),
    FP           = col_double(),
    FN           = col_double(),
    TN           = col_double(),
    Total        = col_double()
  )
)

# -----------------------------------------------------------
# 4) Filter: Nur Neighborhoods, die TP > 0 oder FN > 0 haben
# -----------------------------------------------------------
if (filter_flag) {
  df <- df %>% filter(TP > 0 | FN > 0)
}

# -----------------------------------------------------------
# 5) In Long‐Format transformieren für geom_col (gestapelte Balken)
# -----------------------------------------------------------
df_long <- df %>%
  select(Neighborhood, TP, FP, FN, TN, Total) %>%
  pivot_longer(
    cols      = c(TP, FP, FN, TN),
    names_to  = "ConfusionType",
    values_to = "Count"
  ) %>%
  mutate(Neighborhood = factor(Neighborhood, levels = unique(Neighborhood)))

# Totals separat, damit wir “Total”-Labels über den Balken positionieren:
totals_df <- df %>%
  select(Neighborhood, Total) %>%
  mutate(Neighborhood = factor(Neighborhood, levels = levels(df_long$Neighborhood)))

# -----------------------------------------------------------
# 6) Wir fixieren x‐Achse auf [0 … 5000]
# -----------------------------------------------------------
max_x <- 500
message("[INFO] Fixe x-Achse von 0 bis ", max_x)

# -----------------------------------------------------------
# 7) Physikalische Maße definieren (Balkenhöhe 0,25 cm, Zwischenraum 0,15 cm)
# -----------------------------------------------------------
inch_per_cm <- 1/2.54

# Balkenhöhe = 0,25 cm = 0,098425"
bar_height_in <- 0.25 * inch_per_cm

# Zwischenraum = 0,15 cm = 0,059055"
inter_bar_gap_in <- 0.15 * inch_per_cm

# Kategorie-Höhe = Balkenhöhe + Zwischenraum = 0,25 + 0,15 = 0,40 cm → in Zoll:
category_height_in <- bar_height_in + inter_bar_gap_in  # ≈ 0.15748"

# Bar-Width-Faktor: width × category_height_in = Balkenhöhe.  
# Wir wollen Balkenhöhe = bar_height_in  ⇒ width = bar_height_in / category_height_in
bar_width <- bar_height_in / category_height_in  # ≈ 0.625

n_labels <- length(unique(df_long$Neighborhood))
message("[INFO] Anzahl aller Neighborhood-Zeilen nach Filter: ", n_labels)

# DIN A3-Quer‐Maße in Zoll (420 mm × 297 mm → 16.535" × 11.693")
a3_width_in  <- 420 / 25.4   # ≈ 16.535"
a3_height_in <- 297 / 25.4   # ≈ 11.693"
message(sprintf("[INFO] A3-Quer = %.3f\" × %.3f\" (%.1f × %.1f cm)", 
                a3_width_in, a3_height_in, a3_width_in * 2.54, a3_height_in * 2.54))

# Zeichenbereich vertical: a3_height_in − 1" (oben/unten je 0,5")
draw_area_in <- a3_height_in - 1  # ≈ 10.693"

# Kategorien pro Seite:
n_per_page <- floor(draw_area_in / category_height_in)
message("[INFO] Pro A3-Seite passen bis zu ", n_per_page, " Kategorien (Neighborhoods).")

# Anzahl Seiten:
n_pages <- ceiling(n_labels / n_per_page)
message("[INFO] Benötigte A3-Seiten: ", n_pages)

# -----------------------------------------------------------
# 8) Feste linke Margin: 10 mm
# -----------------------------------------------------------
left_margin_mm <- 10
message("[INFO] Linke Margin fest = ", left_margin_mm, " mm")

# -----------------------------------------------------------
# 9) “Total” für y-Position cappen (alle > max_x → max_x) + Offset
# -----------------------------------------------------------
totals_df <- totals_df %>%
  mutate(
    PlotY = ifelse(Total > max_x, max_x, Total),
    PlotY = PlotY + max_x * 0.01  # 1 % über max_x
  )

# -----------------------------------------------------------
# 10) Funktion: Erzeuge Balken-Plot für eine Teilmenge von Neighborhoods
# -----------------------------------------------------------
make_bar_page <- function(neighborhood_subset, page_index) {
  sub_long   <- df_long   %>% filter(Neighborhood %in% neighborhood_subset)
  sub_totals <- totals_df %>% filter(Neighborhood %in% neighborhood_subset)
  
  sub_long$Neighborhood   <- factor(sub_long$Neighborhood, levels = neighborhood_subset)
  sub_totals$Neighborhood <- factor(sub_totals$Neighborhood, levels = neighborhood_subset)
  
  page_title <- if (n_pages == 1) {
    "Neighborhood Confusion"
  } else {
    sprintf("Neighborhood Confusion (Seite %d von %d)", page_index, n_pages)
  }
  
  p <- ggplot(sub_long, aes(x = Neighborhood, y = Count, fill = ConfusionType)) +
    # 1) Gestapelte Balken ohne Rand, Breite = bar_width
    geom_col(width = bar_width) +
    
    # 2) Horizontaler Flip, fixe x-Skala [0, max_x], und Clip deaktivieren
    coord_flip(ylim = c(0, max_x), clip = "off") +
    
    # 3) Farben für TP/FP/FN/TN
    scale_fill_manual(
      values = c(
        "TP" = "#1b9e77",
        "FP" = "#d95f02",
        "FN" = "#7570b3",
        "TN" = "#e7298a"
      ),
      name = "Confusion Type"
    ) +
    
    # 4) Theme mit 12 pt, dabei rechte Margin auf 10 mm vergrößern
    theme_minimal(base_size = 12) +
    theme(
      # y-Achsentext (Neighborhood) in 12 pt, rechtsbündig
      axis.text.y         = element_text(size = 12, hjust = 1, margin = margin(r = 0, unit = "mm")),
      axis.title.y        = element_blank(),
      # x-Achsentext (Count) in 12 pt
      axis.title.x        = element_text(size = 12),
      axis.text.x         = element_text(size = 12),
      # Legende nur auf Seite 1
      legend.position     = if (page_index == 1) "top" else "none",
      # keine horizontale Gitterlinien zwischen den Balken
      panel.grid.major.y  = element_blank(),
      # Linke Margin fest = left_margin_mm mm; rechte Margin auf 10 mm
      plot.margin         = unit(c(5, 10, 5, left_margin_mm), "mm")
    ) +
    ylab("Absolute Count") +
    xlab(NULL) +
    ggtitle(page_title) +
    
    # 5) “Total”-Beschriftung in 12 pt, Labels werden rechts von der Bar gezeichnet
    geom_text(
      data          = sub_totals,
      inherit.aes   = FALSE,
      aes(x = Neighborhood, y = PlotY, label = Total),
      color         = "black",
      size          = 12 * 0.3527778,  # 1 pt = 0.3527778 mm → 12 pt ≈ 4.233 mm
      hjust         = 0,
      check_overlap = TRUE
    ) +
    
    # 6) 10 % Headroom oberhalb von max_x (statt 5 %)
    expand_limits(y = max_x * 1.10)
  
  return(p)
}

# -----------------------------------------------------------
# 11) Splitte alle Neighborhood-Levels in Blöcke à n_per_page
# -----------------------------------------------------------
all_levels <- levels(df_long$Neighborhood)
pages_list <- split(all_levels, ceiling(seq_along(all_levels) / n_per_page))
# pages_list: Liste mit Länge n_pages; jede Komponente = Vektor von Neighborhood-Levels.

# -----------------------------------------------------------
# 12) Erzeuge Multi-Page A3-PDF (je Seite ein Plot)
# -----------------------------------------------------------
pdf(
  file      = output_file,
  width     = a3_width_in,   # 16.535" (A3-Quer-Breite)
  height    = a3_height_in,  # 11.693" (A3-Quer-Höhe)
  onefile   = TRUE,
  paper     = "special"
)

for (i in seq_along(pages_list)) {
  this_subset <- pages_list[[i]]
  p_page <- make_bar_page(this_subset, i)
  print(p_page)
}

dev.off()

message("✔ Fertig: A3-Quer-PDF geschrieben → ", output_file)
