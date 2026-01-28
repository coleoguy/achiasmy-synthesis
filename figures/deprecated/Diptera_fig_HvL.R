#Diptera Fig
#Blackmon lab
#Meg McConnell (22/12/25)
setwd("~/GitHub/achiasmy-synthesis/data")

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

#Higher vs. Lower Diptera
set.seed(1)  # reproducible random pick among duplicates

# ---- 1) Read data ----
dat <- read_csv("Diptera achiasmy.csv", show_col_types = FALSE)

# ---- 2) Clean taxonomy text ----
tax_cols <- c("Class", "Order", "Suborder", "Family", "Genus", "Species")
dat <- dat %>%
  mutate(across(any_of(tax_cols), ~ str_squish(as.character(.x))))

# ---- 3) Make Haploid + organism key ----
dat <- dat %>%
  mutate(
    `Diploid Number` = suppressWarnings(as.numeric(`Diploid Number`)),
    Haploid = `Diploid Number` / 2,
    organism_id = paste(Class, Order, Suborder, Family, Genus, Species, sep = "|")
  ) %>%
  filter(!is.na(Haploid), !is.na(Suborder), Suborder != "")

# ---- 4) Randomly keep 1 row per unique organism ----
dat_unique <- dat %>%
  group_by(organism_id) %>%
  slice_sample(n = 1) %>%
  ungroup()

# ---- 5) Rename suborders for plotting/legend ----
dat_unique2 <- dat_unique %>%
  mutate(
    Suborder_plot = case_when(
      Suborder == "Brachycera" ~ "Higher Diptera",
      Suborder == "Nematocera" ~ "Lower Diptera",
      TRUE ~ Suborder
    )
  )

# ---- 6) n per suborder (after de-dup & rename) -> legend labels ----
sub_n <- dat_unique2 %>%
  count(Suborder_plot, name = "n") %>%
  arrange(Suborder_plot)

legend_labels <- setNames(
  paste0(sub_n$Suborder_plot, " (n=", sub_n$n, ")"),
  sub_n$Suborder_plot
)

# ---- 7) Frequency per haploid per suborder ----
count_df <- dat_unique2 %>%
  count(Haploid, Suborder_plot, name = "Frequency") %>%
  filter(!is.na(Haploid), !is.na(Suborder_plot))

# ---- 8) Fill missing combos so dodged bars stay same width ----
count_df2 <- count_df %>%
  complete(
    Haploid       = seq(2, max(Haploid, na.rm = TRUE), by = 1),
    Suborder_plot = unique(Suborder_plot),
    fill = list(Frequency = 0)
  )

# ---- 9) Plot ----
p <- ggplot(count_df2, aes(x = Haploid, y = Frequency, fill = Suborder_plot)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.85) +
  scale_x_continuous(breaks = seq(2, max(count_df2$Haploid, na.rm = TRUE), by = 1)) +
  labs(x = "Haploid chromosome number", y = "Frequency", fill = "Suborder") +
  scale_fill_discrete(breaks = names(legend_labels), labels = unname(legend_labels)) +
  theme_classic()

print(p)

# Optional save:
# ggsave("haploid_frequency_by_suborder.png", p, width = 7.5, height = 7, dpi = 300)

