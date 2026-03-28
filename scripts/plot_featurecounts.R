library(tidyverse)

# ---- Read featureCounts summary ----
fc_summary <- read.delim("counts/gene_counts.txt.summary", check.names = FALSE)

# ---- Sample mapping: SRR -> condition labels ----
sample_map <- c(
  "SRR5305524" = "Control_82",
  "SRR5305480" = "Control_81",
  "SRR5305508" = "Control_82",
  "SRR5305483" = "Control_73",
  "SRR5305576" = "Control_97",
  "SRR5305546" = "AD_82",
  "SRR5305567" = "AD_82",
  "SRR5305591" = "AD_81",
  "SRR5305558" = "AD_73",
  "SRR5305477" = "AD_97"
)

# ---- Prep summary for plotting ----
plot_df <- fc_summary %>%
  pivot_longer(-Status, names_to = "sample", values_to = "reads") %>%
  mutate(
    sample = basename(sample),
    sample = sub("\\.Aligned\\.sortedByCoord\\.out\\.bam$", "", sample),
    sample = recode(sample, !!!sample_map)
  ) %>%
  filter(reads > 0) %>%
  mutate(
    group = ifelse(grepl("^AD", sample), "AD", "Control"),
    sample = factor(sample, levels = rev(c(
      "Control_73", "Control_81", "Control_82",
      "Control_82", "Control_97",
      "AD_73", "AD_81", "AD_82",
      "AD_82", "AD_97"
    )))
  )

# ---- Plot ----
ggplot(plot_df, aes(x = sample, y = reads, fill = Status)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "featureCounts Assignment Summary (10 Samples)",
    subtitle = "GSE95587: Alzheimer's Disease vs Control (fusiform gyrus)",
    x = "Samples",
    y = "# Reads"
  ) +
  theme_gray(base_size = 12)

ggsave("counts/featurecounts_summary.pdf", width = 10, height = 6)
ggsave("counts/featurecounts_summary.png", width = 10, height = 6, dpi = 150)
