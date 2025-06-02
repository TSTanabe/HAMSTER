#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(grid)  # Für unit()
})

# ------------------------------------------------------------
# Kommandozeilenargumente einlesen:
#   args[1] = Pfad zur Input-TSV
#   args[2] = Pfad/Name der Output-PDF
#   args[3] = Spaltenname für Y-Achse (z. B. "MCC")
# ------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop(
    "Bitte drei Argumente übergeben:\n",
    "  1) Input-TSV-Datei\n",
    "  2) Output-PDF-Dateiname/-Pfad\n",
    "  3) Spaltenname für Y-Achse (z. B. 'MCC')\n",
    "Beispiel:\n",
    "  Rscript Plot_performance_matrix.R all_performance.txt output.pdf MCC"
  )
}
input_file <- args[1]
output_pdf <- args[2]
y_col      <- args[3]

# ------------------------------------------------------------
# Prüfen, ob Input-TSV existiert
# ------------------------------------------------------------
if (!file.exists(input_file)) {
  stop("Input file '", input_file, "' nicht gefunden.")
}

# ------------------------------------------------------------
# Einlesen der TSV-Datei
#  - erste Spalte: character (Name_prefix)
#  - alle anderen Spalten: automatisch (meist numeric)
# ------------------------------------------------------------
df <- read_tsv(input_file, col_types = cols())

# Name der ersten Spalte (= X-Achse)
x_col <- colnames(df)[1]

# ------------------------------------------------------------
# Prüfen, ob y_col existiert
# ------------------------------------------------------------
if (!(y_col %in% colnames(df))) {
  stop("Spalte '", y_col, "' nicht in der TSV-Datei gefunden.")
}

# ------------------------------------------------------------
# Prefix und Label aus der ersten Spalte extrahieren:
#   prefix = alles vor dem ersten "_"
#   label  = alles nach dem ersten "_"
# ------------------------------------------------------------
df <- df %>%
  mutate(
    raw_name = .data[[x_col]],
    prefix   = sub("^(.*?)_.*$", "\\1", raw_name),
    label    = sub("^[^_]+_(.*)$", "\\1", raw_name)
  )

# ------------------------------------------------------------
# Labels nach maximalem Y-Wert sortieren, damit die x-Achse nach
# Gruppen mit den höchsten Werten angeordnet ist
# ------------------------------------------------------------
label_order <- df %>%
  group_by(label) %>%
  summarize(maxVal = max(.data[[y_col]], na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(maxVal)) %>%
  pull(label)

df <- df %>%
  mutate(label = factor(label, levels = label_order))

# ------------------------------------------------------------
# Dynamische PDF-Breite basierend auf der Anzahl der Labels:
#   - Basis-Breite: 8 Zoll
#   - Pro Label: 0.3 Zoll
# ------------------------------------------------------------
n_labels        <- nlevels(df$label)
width_per_label <- 0.3     # Zoll pro Label
base_width      <- 8       # Mindestbreite in Zoll
pdf_width       <- max(base_width, n_labels * width_per_label)
pdf_height      <- 5       # feste Höhe

# ------------------------------------------------------------
# Falls Ausgabedatei schon existiert, Abbruch (nicht überschreiben)
# ------------------------------------------------------------
if (file.exists(output_pdf)) {
  message("[INFO] Output file '", output_pdf, "' existiert bereits – Skript wird beendet.")
  quit(status = 0)
}

# ------------------------------------------------------------
# PDF öffnen
# ------------------------------------------------------------
pdf(output_pdf, width = pdf_width, height = pdf_height)

# ------------------------------------------------------------
# Plot erstellen:
#   - X-Achse: label (ohne Prefix)
#   - Y-Achse: y_col
#   - Farbe/Form: prefix
#   - Hilfslinien auf beiden Achsen (major grid lines)
#   - X-Achse: Labels schräg, nicht überfüllt dank guide_axis
# ------------------------------------------------------------
ggplot(df, aes(x = label, y = .data[[y_col]], color = prefix, shape = prefix)) +
  geom_point(size = 2) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.x = element_line(color = "grey80", size = 0.5),
    panel.grid.major.y = element_line(color = "grey80", size = 0.5),
    panel.grid.minor   = element_blank(),
    # X-Achsentext schräg und kleiner
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 8, color = "black"),
    axis.text.y      = element_text(size = 8, color = "black"),
    # Achsentitel kleiner
    axis.title.x     = element_text(size = 9, color = "black"),
    axis.title.y     = element_text(size = 9, color = "black"),
    # Legende kleiner
    legend.title     = element_text(size = 9, color = "black"),
    legend.text      = element_text(size = 8, color = "black"),
    # Plot-Titel
    plot.title       = element_text(size = 11, color = "black", hjust = 0.5),
    plot.margin      = unit(c(5, 5, 5, 5), "mm")
  ) +
  # X-Achse mit guide_axis, um Tickmarks nicht zu überfüllen
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE, n.dodge = 1)) +
  xlab("HMM") +
  ylab(y_col) +
  ggtitle("self-recognition")

dev.off()

message("✔ PDF gespeichert unter: ", output_pdf,
        " (Breite=", round(pdf_width, 2), "\" × Höhe=", pdf_height, "\")")
