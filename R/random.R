library(ggplot2)
library(dplyr)

df <- data.frame(
  pathway = c(
    "Metabolic pathways",
    "ABC transporters",
    "Lipopolysaccharide biosynthesis",
    "Biosynthesis of unsaturated fatty acids",
    "Cationic antimicrobial peptide (CAMP) resistance",
    "Bacterial secretion system",
    "Biosynthesis of cofactors",
    "Thiamine metabolism",
    "Protein export"
  ),
  count = c(2, 2, 1, 1, 1, 1, 1, 1, 1)
)

df <- df %>% arrange(count)

bar <- ggplot(df, aes(x = count, y = reorder(pathway, count), fill = count)) +
  geom_col() +
  geom_text(aes(label = count), hjust = -0.15, size = 4) +
  scale_fill_gradient(low = "#1B98E026", high="#225ea8") +
  labs(
    x = "Number of associated enriched KOs",
    y = NULL,
    fill = "Count"
  ) +
  theme_classic(base_size = 12) +
  xlim(0, max(df$count) + 0.5)

ggsave(
  "results/aim3/bar.png",
  plot = bar,
  width = 10, height = 5, units = "in", dpi = 300
)
bar