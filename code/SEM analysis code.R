#Script for SEM analysis 

library(lme4)
library(lmerTest)
library(emmeans)
library(ggbeeswarm)
library(patchwork)
library(tidyverse)
library(janitor)
library(stringr)
library(readr)

files <- list.files(
  pattern = "SEM_analysis\\.csv$",
  full.names = TRUE
)

length(files) 


raw <- files %>%
  set_names(basename(.)) %>%
  purrr::map_dfr(~ readr::read_csv(.x, show_col_types = FALSE),
                 .id = "filename")

data <- raw %>%
  select(-`...1`) %>%
  janitor::clean_names() %>%
  rename(feret_um = feret) %>%
  mutate(
    replicate = stringr::str_extract(filename, "Sample[A-B]"),
    image_id  = stringr::str_match(filename, "Sample[A-B]_([0-9]+(?:\\.[0-9]+)?)")[,2],
    roi_id    = stringr::str_match(filename, "\\.([abc])_SEM_particlesFIJI\\.csv$")[,2],
    roi_id    = if_else(is.na(roi_id), "single", roi_id),
    replicate = factor(replicate),
    image_id  = factor(image_id),
    roi_id    = factor(roi_id, levels = c("a","b","c","single"))
  )

data <- data %>%
  mutate(
    image_id = if_else(filename == "SampleA_0.11.b_SEM_particlesFIJI.csv", "11", as.character(image_id)),
    roi_id   = if_else(filename == "SampleA_0.11.b_SEM_particlesFIJI.csv", "b",  as.character(roi_id)),
    image_id = factor(image_id),
    roi_id   = factor(roi_id, levels = c("a","b","c","single"))
  )

count(data, replicate, image_id, roi_id) %>%
  arrange(replicate, suppressWarnings(as.numeric(as.character(image_id))), roi_id) %>%
  print(n = 200)

roi_summary <- data %>%
  group_by(replicate, image_id, roi_id, filename) %>%
  summarise(
    n_objects    = n(),
    mean_feret   = mean(feret_um, na.rm = TRUE),
    median_feret = median(feret_um, na.rm = TRUE),
    p90_feret    = quantile(feret_um, 0.9, na.rm = TRUE),
    max_feret    = max(feret_um, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(low_count = n_objects < 5)

roi_summary %>% arrange(replicate, as.numeric(as.character(image_id)), roi_id) %>% print(n = 200)

roi_summary_5plus <- roi_summary %>% filter(n_objects >= 5)
nrow(roi_summary); nrow(roi_summary_5plus)

roi_summary %>%
  filter(n_objects < 5) %>%
  arrange(replicate, as.numeric(as.character(image_id)), roi_id)
view(roi_summary)

p_p90 <- ggplot(roi_summary, aes(replicate, p90_feret)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(
    width = 0.15,
    size = 2,
    aes(shape = low_count)
  ) +
  labs(
    x = "Migration tube (biological replicate)",
    y = "90th percentile Feret diameter (µm)"
  ) +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    legend.position = "none",
    text = element_text(family = "Arial"),
    axis.title.x = element_text(size = 17),
    axis.title.y = element_text(size = 17),
    axis.text.x  = element_text(size = 15),
    axis.text.y  = element_text(size = 15)
  )


p_median<-ggplot(roi_summary, aes(replicate, median_feret)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(
    width = 0.15,
    size = 2,
    aes(shape = low_count)
  ) +
  labs(
    x = "Migration tube",
    y = "Median Feret diameter (µm)",
    shape = "Low object count"
  ) +
  theme_classic()


image_summary <- roi_summary %>%
  group_by(replicate, image_id) %>%
  summarise(
    p90_feret_img = median(p90_feret),
    median_feret_img = median(median_feret),
    .groups = "drop"
  )

p_image<-ggplot(image_summary, aes(replicate, p90_feret_img)) +
  geom_point(size = 3) +
  labs(
    x = "Migration tube",
    y = "Image-level p90 Feret diameter (µm)"
  ) +
  theme_classic()



ggsave(
  filename = "Figure_p90_Feret_ROI.jpeg",
  plot     = p_p90,
  device   = "jpeg",
  width    = 6,
  height   = 5,
  units    = "in",
  dpi      = 900
)

ggsave(
  filename = "Figure_median_Feret_ROI.jpeg",
  plot     = p_median,
  device   = "jpeg",
  width    = 6,
  height   = 5,
  units    = "in",
  dpi      = 900
)

ggsave(
  filename = "Figure_image_level_p90.jpeg",
  plot     = p_image,
  device   = "jpeg",
  width    = 6,
  height   = 5,
  units    = "in",
  dpi      = 900
)
