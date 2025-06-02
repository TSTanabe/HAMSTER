#!/usr/bin/env Rscript

# plot_neighborhood_roc_A4_neu_blackaxes.R
#
# Dieses Skript liest eine TSV-Datei mit mindestens den Spalten „FPR“, „TPR“
# und einer dritten Spalte (Score) ein und erzeugt daraus einen klassischen
# ROC-Plot. Achsenlinien, Achsenbeschriftungen und Hilfslinien werden in Schwarz
# dargestellt. Das Ergebnis wird als einseitiges PDF im DIN A4-Format (Portrait)
# im selben Verzeichnis wie die Eingabedatei abgelegt.
#
# Usage:
#   Rscript plot_neighborhood_roc_A4_neu_blackaxes.R /pfad/zur/datei.tsv
#
# -----------------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(grid)
})

# -----------------------------------------------------------
# 1) Pfad zur TSV-Datei festlegen (hier direkt hartkodiert)
# -----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop(
    "Bitte den Pfad zur TSV-Datei als Argument angeben.\n",
    "Beispiel:\n",
    "  Rscript plot_neighborhood_roc_A4_labels_final.R /pfad/zur/datei.tsv"
  )
}
tsv_path <- args[1]
# -----------------------------------------------------------
# 2) Prüfen, ob die Datei existiert
# -----------------------------------------------------------
if (!file.exists(tsv_path)) {
  stop("Datei nicht gefunden: ", tsv_path)
}

# -----------------------------------------------------------
# 3) Ausgabepfad im gleichen Verzeichnis bestimmen
# -----------------------------------------------------------
input_dir   <- dirname(tsv_path)
input_base  <- tools::file_path_sans_ext(basename(tsv_path))
output_file <- file.path(input_dir, paste0(input_base, "_roc.pdf"))

# Wenn das PDF schon existiert, Skript beenden
if (file.exists(output_file)) {
  message("[INFO] Ziel-PDF existiert bereits: ", output_file, "\n",
          "[INFO] Abbruch, keine neue Datei angelegt.")
  quit(status = 0)
}
message("[INFO] Ausgabe-PDF: ", output_file)

# -----------------------------------------------------------
# 4) TSV-Datei einlesen (FPR, TPR und Score erwarten)
# -----------------------------------------------------------
message("[INFO] Lese TSV-Datei ein: ", tsv_path)
df_raw <- tryCatch(
  read_tsv(
    tsv_path,
    col_types = cols(
      FPR    = col_double(),
      TPR    = col_double(),
      .default = col_double()
    )
  ),
  error = function(e) {
    stop("Fehler beim Einlesen der TSV-Datei: ", e$message)
  }
)

# -----------------------------------------------------------
# 5) Prüfen, ob Spalten FPR und TPR vorhanden sind
# -----------------------------------------------------------
if (!all(c("FPR", "TPR") %in% names(df_raw))) {
  stop("Die TSV-Datei muss die Spalten »FPR« und »TPR« enthalten.")
}

# -----------------------------------------------------------
# 6) Score-Spalte bestimmen (Dritte Spalte)
# -----------------------------------------------------------
if (ncol(df_raw) < 3) {
  stop("Die TSV-Datei muss mindestens drei Spalten haben (FPR, TPR und Score).")
}
score_colname <- names(df_raw)[3]
message("[INFO] Verwende Spalte '", score_colname, "' als Score.")

# -----------------------------------------------------------
# 7) Daten für den ROC-Plot aufbereiten
# -----------------------------------------------------------
df <- df_raw %>%
  select(FPR, TPR, Score = !!sym(score_colname)) %>%
  # Nur Zeilen mit gültigen Werten behalten
  filter(
    !is.na(FPR), !is.na(TPR), !is.na(Score),
    FPR >= 0, FPR <= 1,
    TPR >= 0, TPR <= 1
  ) %>%
  arrange(FPR)

if (nrow(df) == 0) {
  stop("Nach Filterung sind keine gültigen Datenzeilen mehr übrig.")
}
message("[INFO] Anzahl gültiger Zeilen für ROC: ", nrow(df))

# -----------------------------------------------------------
# 8) DIN A4-Abmessungen (Portrait) in Zoll berechnen
# -----------------------------------------------------------
inch_per_cm <- 1/2.54
a4_width_in  <- 210 / 25.4   # ≈ 8.268 Zoll
a4_height_in <- 297 / 25.4   # ≈ 11.693 Zoll
message(sprintf("[INFO] DIN A4 (Portrait) = %.3f\" × %.3f\" (%.1f × %.1f cm)",
                a4_width_in, a4_height_in,
                a4_width_in * 2.54, a4_height_in * 2.54))

# Linke Margin in mm
left_margin_mm <- 10

# -----------------------------------------------------------
# 9) Funktion: ROC-Plot erzeugen (schwarze Achsen, schwarze Hilfslinien)
# -----------------------------------------------------------
make_roc_plot <- function(data) {
  ggplot(data, aes(x = FPR, y = TPR)) +
    # 45°-Diagonale als Referenz (Zufalls-Klassifikator)
    geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
    # ROC-Linie verbinden
    geom_line(color = "#1b9e77", size = 0.8) +
    # Punkte, eingefärbt nach Score
    geom_point(aes(color = Score), size = 2) +
    scale_color_gradient(low = "#d95f02", high = "#7570b3", name = "Score") +
    # Achsen von 0 bis 1
    scale_x_continuous(
      name   = "False Positive Rate (FPR)",
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2)
    ) +
    scale_y_continuous(
      name   = "True Positive Rate (TPR)",
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2)
    ) +
    ggtitle("ROC Curve") +
    theme_minimal(base_size = 12) +
    theme(
      # Achsenlinien und Ticks in Schwarz
      axis.line       = element_line(color = "black"),
      axis.ticks      = element_line(color = "black"),
      # Achsentitel in Schwarz, 12 pt
      axis.title.x    = element_text(color = "black", size = 12),
      axis.title.y    = element_text(color = "black", size = 12),
      # Achsentext (Zahlen) in Schwarz, 12 pt
      axis.text.x     = element_text(color = "black", size = 12),
      axis.text.y     = element_text(color = "black", size = 12),
      # Hilfslinien (major/minor) in Schwarz, gestrichelt
      panel.grid.major = element_line(color = "black", linetype = "dotted"),
      panel.grid.minor = element_line(color = "black", linetype = "dotted"),
      # Hintergrund weiß
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      # Plot-Margin: oben 5 mm, rechts 10 mm, unten 5 mm, links = left_margin_mm
      plot.margin     = unit(c(5, 10, 5, left_margin_mm), "mm"),
      # Legende rechts
      legend.position = "right"
    )
}

# -----------------------------------------------------------
# 10) Plot erzeugen und als PDF speichern
# -----------------------------------------------------------
message("[INFO] Erzeuge ROC-Plot …")
roc_plot <- make_roc_plot(df)

message("[INFO] Speichere PDF …")
ggsave(
  filename = output_file,
  plot     = roc_plot,
  width    = a4_width_in,
  height   = a4_height_in,
  units    = "in",
  device   = "pdf"
)

message("✔ Fertig: ROC-PDF erstellt → ", output_file)
