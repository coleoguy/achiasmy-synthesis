# Mammal Fig (Subclass split)
# Blackmon lab
# Meg McConnell (22/12/25)

setwd("~/GitHub/achiasmy-synthesis/data")

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

set.seed(1)  # reproducible random choices

# ---- 1) Read data ----
dat <- read_csv("Mammalia achiasmy.csv", show_col_types = FALSE)

# ---- 2) Clean taxonomy text ----
tax_cols <- c("Class", "Subclass", "Order", "Family", "Genus", "Species")
dat <- dat %>%
  mutate(across(any_of(tax_cols), ~ str_squish(as.character(.x))))

# ---- 3) Make Haploid + organism key ----
dat <- dat %>%
  mutate(
    `Diploid Number` = suppressWarnings(as.numeric(`Diploid Number`)),
    Haploid     = `Diploid Number` / 2,
    organism_id = paste(Class, Subclass, Order, Family, Genus, Species, sep = "|")
  ) %>%
  filter(!is.na(Haploid), !is.na(Subclass), Subclass != "")

# ---- 4) Randomly keep 1 row per unique organism (duplicates) ----
dat_unique <- dat %>%
  group_by(organism_id) %>%
  slice_sample(n = 1) %>%
  ungroup()

# ---- 5) Rename subclasses for plotting/legend ----
dat_unique2 <- dat_unique %>%
  mutate(
    Subclass_plot = case_when(
      Subclass == "Metatheria" ~ "Marsupials",
      Subclass == "Theria"     ~ "Placental Mammals",
      TRUE ~ Subclass
    )
  )

# ---- 6) n per subclass (after de-dup & rename) -> legend labels ----
sub_n <- dat_unique2 %>%
  count(Subclass_plot, name = "n") %>%
  arrange(Subclass_plot)

legend_labels <- setNames(
  paste0(sub_n$Subclass_plot, " (n=", sub_n$n, ")"),
  sub_n$Subclass_plot
)

# ---- 7) Convert haploid to INTEGER bins:
#        - if Haploid is already integer -> keep it
#        - if Haploid is x.5 (e.g., 19.5) -> randomly pick floor (19) OR floor-1 (18)
#          (i.e., 0.5 rounds DOWN either 0 or 1 extra)
dat_unique2 <- dat_unique2 %>%
  mutate(
    Haploid_int = case_when(
      abs(Haploid - round(Haploid)) < 1e-8 ~ as.integer(round(Haploid)),
      abs(Haploid*2 - round(Haploid*2)) < 1e-8 ~ {  # halves like 19.5
        f <- floor(Haploid)                        # 19
        # randomly choose f or f-1
        as.integer(f - sample(c(0L, 1L), size = 1))
      },
      TRUE ~ NA_integer_
    )
  ) %>%
  filter(!is.na(Haploid_int), Haploid_int >= 2)

# ---- 8) Frequency per integer haploid per subclass ----
sub_levels <- c("Marsupials", "Placental Mammals")
count_df <- dat_unique2 %>%
  mutate(Subclass_plot = factor(Subclass_plot, levels = sub_levels)) %>%
  count(Haploid_int, Subclass_plot, name = "Frequency")

# ---- 9) Fill missing combos so dodged bars stay same width ----
maxH <- max(count_df$Haploid_int, na.rm = TRUE)

count_df2 <- count_df %>%
  complete(
    Haploid_int   = seq(2, maxH, by = 1),
    Subclass_plot = factor(sub_levels, levels = sub_levels),
    fill = list(Frequency = 0)
  ) %>%
  group_by(Haploid_int, Subclass_plot) %>%                 # ensure 1 row per combo
  summarise(Frequency = sum(Frequency), .groups = "drop") %>%
  mutate(Haploid_f = factor(Haploid_int, levels = seq(2, maxH, by = 1)))

# ---- 10) Plot ----
p <- ggplot(count_df2, aes(x = Haploid_f, y = Frequency, fill = Subclass_plot)) +
  geom_col(
    position = position_dodge2(width = 0.9, preserve = "single"),
    width = 0.85
  ) +
  labs(
    x = "Haploid chromosome number",
    y = "Frequency",
    fill = "Subclass"
  ) +
  scale_fill_discrete(
    breaks = names(legend_labels),
    labels = unname(legend_labels)
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

print(p)

# Optional save:
# ggsave("haploid_frequency_by_subclass.png", p, width = 9, height = 6, dpi = 300)
