#' Supplementary code to the article: 
#' Klinkovska et al. (2025) Half a century of temperate non-forest vegetation changes: 
#' no net loss in species richness, but considerable shifts in taxonomic and functional composition. Global Change Biology.
#' 
#' Author: Klara Klinkovska, Helge Bruelheide 2024-09-24
#' R version 4.3.2

library(broom) # version 1.0.5
library(mgcv) # version 1.9-1
library(vegan) # version 2.6-8
library(ggrepel) # version 0.9.4
library(patchwork) # version 1.2.0
library(tidyverse) # version 2.0.0
# dplyr version 1.1.3
# forcats version 1.0.0
# ggplot2 version 3.5.1
# purrr version 1.0.2
# readr version 2.1.4
# stringr version 1.5.1
# tibble version 3.2.1
# tidyr version 1.3.1

# color settings
veg_col <- c('A' = "#F8766D", 'M' = "#D39200", 'R' = "#DB72FB", 'TD_wet' = "#619CFF", 'TD_mesic' = "#00B9E3", 
             'TE' = "#00C19F", 'TFG' = "#00BA38", 'TH' = "#A3A500", 'X' = "#FF61CC")

# load data ---------------------------------------------------------------

# original data
head <- read_csv("data/Klinkovska_et_al_half_a_century_of_temperate_vegetation_change_head.csv") 

# different plot sizes
diff_area <-  read_csv('data/diff_area.csv') |> 
  select(RS_CODE_PLOT = RS_CODE_PLOT.x, from.PlotObservationID, to.PlotObservationID, 
         Releve_area_m2.x, Releve_area_m2.y) |> 
  mutate(prop_diff = Releve_area_m2.x/Releve_area_m2.y)

mean(diff_area$prop_diff)
min(diff_area$prop_diff)
max(diff_area$prop_diff)


# interval changes 
interval_change <- read_csv('results/intervals_all.csv') |> 
  mutate(ts_length = to-from) |> 
  mutate(prop_diff = to_area/from_area, # exclude plots of different sizes for diversity indices
         interval_change = if_else(variable %in% c('richness', 'shannon') & prop_diff != 1, 
                                  NA, interval_change, interval_change), 
         interval_change = if_else(variable == 'evenness' & (prop_diff >= 2 | prop_diff <= 0.5), 
                                  NA, interval_change, interval_change), 
         interval_change_sign = sign(interval_change), 
         test = 'interval') |> 
  select(RS_CODE_PLOT, from.PlotObservationID, to.PlotObservationID, from, to, ts_length, 
         variable, veg_type0, interval_change, test, Releve_area_m2, from_area, to_area, prop_diff,
         interval_change_sign, RS_CODE, Longitude, Latitude) 

plot_char <- read_csv('results/plot_characteristics.csv')

# changes estimated with linear models 
plot_change_lm <- read_csv('results/lm_trends.csv') |> 
  rename('lm_change' = 'est') |> 
  left_join(diff_area |> distinct(RS_CODE_PLOT, .keep_all = T) |> # to exclude species richness and shannon index changes for observations of different sizes
              select(RS_CODE_PLOT, prop_diff), 
            by = 'RS_CODE_PLOT') |> 
  left_join(diff_area |> distinct(RS_CODE_PLOT, .keep_all = T) |> # to exclude evenness changes for observations of different sizes
              filter(prop_diff >= 2 | prop_diff <= 0.5) |> 
              select(RS_CODE_PLOT, prop_diff), 
            by = 'RS_CODE_PLOT') |> 
  mutate(lm_change = ifelse(variable %in% c('richness', 'shannon') & !is.na(prop_diff.x), NA, lm_change), 
         lm_change = ifelse(variable %in% c('evenness') & !is.na(prop_diff.y), NA, lm_change),
         lm_change_sign = sign(lm_change)) |> 
  left_join(interval_change |> 
              filter(test == 'interval' & variable == 'richness') |> 
              group_by(RS_CODE_PLOT) |> 
              summarize(n_surveys = n() + 1) |>
              ungroup()) |> 
  select(-prop_diff.x, -prop_diff.y) |> 
  left_join(head |> group_by(RS_CODE_PLOT) |>
              mutate(min_year = min(year),
                     max_year = max(year),
                     diff_year = max_year - min_year,
                     mid_year = min_year + diff_year/2) |>
              ungroup() |> 
              distinct(RS_CODE_PLOT, .keep_all = T)) |> 
  mutate(RS_CODE = as.factor(RS_CODE))

# order traits
trait_order <- read_csv('data/variable_order.csv')

# gam binomial model ----------------------------------------------------------
# linear change approach by habitat
plot_change_lm2 <- plot_change_lm |>
  filter(lm_change_sign != 0 & str_detect(variable, 'CR_', negate = T)) |> 
  mutate(lm_change_sign = ifelse(lm_change_sign == -1, 0, lm_change_sign)) |> 
  nest(data = -variable) |> 
  mutate(m_tidy = map(data, ~gam(lm_change_sign ~ veg_type0 - 1 + s(Latitude, Longitude, bs = "sos") + s(RS_CODE, bs = "re"), 
                          weights = sqrt(n_surveys),  family = "binomial", method = "REML", data = .x) |> 
                        tidy(conf.int = T, parametric = T)), 
         n = map(data, ~.x |> group_by(veg_type0) |> count())) |> 
  unnest(c(m_tidy, n)) |> 
  mutate(across(c(estimate, conf.low, conf.high), ~1/(1+exp(-.)))) |> 
  select(-data)

write_csv(plot_change_lm2, 'results/gam_lm_change_binomial.csv')

# linear change approach, all habitats together
plot_change_lm_all <- plot_change_lm |>
  filter(lm_change_sign != 0 & str_detect(variable, 'CR_', negate = T)
         ) |> 
  mutate(lm_change_sign = ifelse(lm_change_sign == -1, 0, lm_change_sign)) |> 
  nest(data = -variable) |> 
  mutate(m_tidy = map(data, ~gam(lm_change_sign ~ 1 + s(Latitude, Longitude, bs = "sos") + s(RS_CODE, bs = "re"), 
                                 weights = sqrt(n_surveys),  family = "binomial", method = "REML", data = .x) |> 
                        tidy(conf.int = T, parametric = T)), 
         n = map(data, ~.x |> 
                   count())) |> 
  unnest(c(m_tidy, n)) |> 
  mutate(across(c(estimate, conf.low, conf.high), ~1/(1+exp(-.)))) |> 
  select(-data)

write_csv(plot_change_lm_all, 'results/gam_lm_change_binomial_all.csv')

# interval change approach by habitat
plot_change_interval2 <- interval_change |>
  filter(interval_change_sign != 0 & test == 'interval' & str_detect(variable, 'CR|shannon|selfing', negate = T)
         ) |> 
  mutate(interval_change_sign = ifelse(interval_change_sign == -1, 0, interval_change_sign), 
         RS_CODE = as.factor(RS_CODE), RS_CODE_PLOT = as.factor(RS_CODE_PLOT)) |> 
  nest(data = -c(variable)) |> 
  mutate(m_tidy = map(data, ~gam(interval_change_sign ~ veg_type0 - 1 + 
                                   s(Latitude, Longitude, bs = "sos") + s(RS_CODE, bs = "re"), 
                          family = "binomial", method = "REML", data = .x) |> 
                        tidy(conf.int = T, parametric = T)), 
         n = map(data, ~.x |> group_by(veg_type0) |> count())) |> 
  unnest(c(m_tidy, n)) |> 
  mutate(across(c(estimate, conf.low, conf.high), ~1/(1+exp(-.)))) |> 
  select(-data)

write_csv(plot_change_interval2, 'results/gam_interval_change_binomial.csv')

# interval change approach, all habitats
plot_change_interval_all <- interval_change |>
  filter(interval_change_sign != 0 & test == 'interval' & 
           str_detect(variable, 'CR_|shannon|selfing', negate = T)) |> 
  mutate(interval_change_sign = ifelse(interval_change_sign == -1, 0, interval_change_sign), 
         RS_CODE = as.factor(RS_CODE), RS_CODE_PLOT = as.factor(RS_CODE_PLOT)) |> 
  nest(data = -c(variable)) |> 
  mutate(m_tidy = map(data, ~gam(interval_change_sign ~ 1 + s(Latitude, Longitude, bs = "sos") + 
                                   s(RS_CODE, bs = "re"), 
                                 family = "binomial", method = "REML", data = .x) |> 
                        tidy(conf.int = T, parametric = T)), 
         n = map(data, ~.x |> 
                   count())) |> 
  unnest(c(m_tidy, n)) |> 
  mutate(across(c(estimate, conf.low, conf.high), ~1/(1+exp(-.)))) |> 
  select(-data)

write_csv(plot_change_interval_all, 'results/gam_interval_change_binomial_all.csv')

# join together 
plot_change_lm2 <- read_csv('results/gam_lm_change_binomial.csv') |> 
  bind_rows(read_csv('results/gam_lm_change_binomial_all.csv'))

plot_change_interval2 <- read_csv('results/gam_interval_change_binomial.csv') |> 
  bind_rows(read_csv('results/gam_interval_change_binomial_all.csv'))

plot_change_all3 <- bind_rows(plot_change_lm2 |> mutate(test = 'lm'), 
                              plot_change_interval2 |> mutate(test = 'interval')) |> 
  mutate(term = str_remove(term, 'veg_type0') |> 
           fct_relevel('(Intercept)', 'A', 'M', 'R', 'TD_wet', 'TD_mesic', 'TE', 'TFG', 'TH', 'X'), 
         test = fct_relevel(test, 'interval', 'lm')) |> 
  filter(variable != 'invsimps')

# gam results forest plot taxonomic diversity ---------------
plot_change_all3 |> 
  filter(str_detect(variable, 'richness|evenness') & test == 'lm'
         ) |>
  mutate(variable = fct_relevel(variable, 'richness', 'evenness')) |> 
  ggplot(aes(x = estimate, xmin = conf.low, xmax = conf.high, y = term, 
             color = term#, 
             #shape = test
             ))+
  geom_vline(xintercept = 0.5, linetype = 'dashed')+
  geom_pointrange(position = position_dodge(width = 0.7))+
  scale_y_discrete(limits = rev,
                   labels = c("Ruderal and \nweed vegetation", 
                              'Dry grasslands',
                              'Sand and shallow-soil \nvegetation',
                              "Nardus grasslands \nand heathlands",
                              'Mesic meadows \nand pastures',
                              "Wet meadows",
                              "Spring and mire \nvegetation",
                              "Wetland vegetation",
                              "Alpine and \nsubalpine vegetation",
                              'All habitats'))+
  scale_color_manual(values = veg_col, guide = 'none')+
  #scale_shape_manual(values = c(17, 16), labels = c('Interval change approach', 'Linear trend approach'))+
  labs(x = 'Probability to increase')+
  facet_wrap(~variable,
             labeller = labeller(variable = c(richness = '(a) Species richness',
                                              evenness = '(b) Pielou evenness')))+
  theme_bw()+
  theme(legend.title = element_blank(), legend.position = 'bottom', strip.background = element_blank(), 
        axis.title.y = element_blank())

ggsave('plots/diversity_lm_binomial.png', width = 7, height = 4)

### summary table
plot_change_sum <- plot_change_all3 |> 
  filter(str_detect(variable, 'wm') & test == 'lm') |> 
  mutate(variable = str_replace(variable, '_wm', '.wm') |> str_replace('_uwm', '.uwm')) |> 
  separate_wider_delim(cols = variable, delim = '.', names = c('variable', 'weighted')) |> 
  left_join(trait_order) |>
  mutate(label = fct_reorder(label, -order), 
         group = fct_reorder(group, order),
         trend = ifelse(p.value < 0.05 & estimate > 0.5, 'increasing', 
                        ifelse(p.value < 0.05 & estimate < 0.5, 'decreasing', NA)), 
         test2 = str_c(test, weighted, sep = '_')) |> 
  filter(!is.na(order))

# write table to appendix
plot_change_all3 |> 
  mutate(variable = str_replace(variable, '_wm', '.wm') |> str_replace('_uwm', '.uwm')) |> 
  separate_wider_delim(cols = variable, delim = '.', names = c('variable', 'weighted'), 
                       too_few = "align_start") |> 
  left_join(trait_order) |>
  mutate(label = fct_reorder(label, -order), 
         group = fct_reorder(group, order), 
         trend = ifelse(estimate < 0.5, 'decreasing', 'increasing')) |> 
  filter(!is.na(order)) |> 
  arrange(order, term, weighted, test) |> 
  select(`Variable category` = group, Variable = label, `Vegetation type` = term, 
         Test = test, Weighted = weighted, Trend = trend, Estimate = estimate, 
         `Standard error` = std.error, `Test statistic` = statistic, 
         `P-value` = p.value, `Lower CI` = conf.low, `Upper CI` = conf.high, 
         `Number of change observations` = n) |> 
  mutate(`Vegetation type` = case_when(`Vegetation type` == '(Intercept)' ~ 'All habitats',
                              `Vegetation type` == 'A' ~ 'Alpine and subalpine vegetation', 
                              `Vegetation type` == 'M' ~ 'Wetland vegetation',
                              `Vegetation type` == 'R' ~ 'Spring and mire vegetation',
                              `Vegetation type` == 'TD_wet' ~ 'Wet meadows',
                              `Vegetation type` == 'TD_mesic' ~ 'Mesic meadows and pastures',
                              `Vegetation type` == 'TE' ~ 'Nardus grasslands and heathlands',
                              `Vegetation type` == 'TFG' ~ 'Sand and shallow-soil vegetation',
                              `Vegetation type` == 'TH' ~ 'Dry grasslands',
                              `Vegetation type` == 'X' ~ 'Ruderal and weed vegetation'), 
         Weighted = ifelse(Weighted == 'uwm', 'unweighted', 'weighted')) |> 
  write_csv('results/gam_interval_lm_binomial_all.csv')

# summary plot
plot_change_sum |> 
  filter(!is.na(trend)) |> 
  ggplot(aes(test2, label, fill = trend, shape = trend))+
  geom_point(size = 2.5)+
  scale_shape_manual(values = c(25, 24))+
  facet_grid(rows = vars(group), cols = vars(term), scales = 'free', space = 'free', switch = 'y',
             labeller = labeller(term = c(A = "Alpine and \nsubalpine vegetation",
                                          M = "Wetland vegetation",
                                          R = "Spring and mire \nvegetation",
                                          TD_mesic = "Mesic meadows \nand pastures",
                                          TD_wet = 'Wet meadows',
                                          TE = "Nardus grasslands \nand heathlands",
                                          TFG = 'Sand and shallow-soil \nvegetation',
                                          TH = 'Dry grasslands',
                                          X = "Ruderal and \nweed vegetation", 
                                          `(Intercept)` = 'All habitats'), 
                                 group = label_wrap_gen(width = 15)))+
  scale_x_discrete(labels = c('unweighted', 'weighted'))+
  theme_bw()+
  theme(legend.title = element_blank(), legend.position = 'bottom', 
        axis.title = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1), 
        strip.text.x = element_text(angle = 90), 
        strip.text.y.left = element_text(angle = 0, face = 'bold'),
        strip.background = element_blank(), 
        panel.spacing.y = unit(0, 'lines'), strip.placement = 'outside')

ggsave('plots/summary_table.png', width = 6.5, height = 13)


### summary table interval version
plot_change_sum <- plot_change_all3 |> 
  filter(str_detect(variable, 'wm') & test == 'interval') |> 
  mutate(variable = str_replace(variable, '_wm', '.wm') |> str_replace('_uwm', '.uwm')) |> 
  separate_wider_delim(cols = variable, delim = '.', names = c('variable', 'weighted')) |> 
  left_join(trait_order) |>
  mutate(label = fct_reorder(label, -order), 
         group = fct_reorder(group, order),
         trend = ifelse(p.value < 0.05 & estimate > 0.5, 'increasing', 
                        ifelse(p.value < 0.05 & estimate < 0.5, 'decreasing', NA)), 
         test2 = str_c(test, weighted, sep = '_')) |> 
  filter(!is.na(order))

# summary plot
plot_change_sum |> 
  filter(!is.na(trend)) |> 
  ggplot(aes(test2, label, fill = trend, shape = trend))+
  geom_point(size = 2.5)+
  scale_shape_manual(values = c(25, 24))+
  facet_grid(rows = vars(group), cols = vars(term), scales = 'free', space = 'free', switch = 'y',
             labeller = labeller(term = c(A = "Alpine and \nsubalpine vegetation",
                                          M = "Wetland vegetation",
                                          R = "Spring and mire \nvegetation",
                                          TD_mesic = "Mesic meadows \nand pastures",
                                          TD_wet = 'Wet meadows',
                                          TE = "Nardus grasslands \nand heathlands",
                                          TFG = 'Sand and shallow-soil \nvegetation',
                                          TH = 'Dry grasslands',
                                          X = "Ruderal and \nweed vegetation", 
                                          `(Intercept)` = 'All habitats'), 
                                 group = label_wrap_gen(width = 15)))+
  scale_x_discrete(labels = c('unweighted', 'weighted'))+
  theme_bw()+
  theme(legend.title = element_blank(), legend.position = 'bottom', 
        axis.title = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1), 
        strip.text.x = element_text(angle = 90), 
        strip.text.y.left = element_text(angle = 0, face = 'bold'),
        strip.background = element_blank(), 
        panel.spacing.y = unit(0, 'lines'), strip.placement = 'outside')

ggsave('plots/summary_table_int.png', width = 6.5, height = 13)

# pca ---------------------------------------------------------------------
plot_change_all4 <- plot_change_all3 |> 
  filter(str_detect(variable, '_uwm') & test == 'lm' & term != '(Intercept)') |>
  mutate(variable = str_remove(variable, '_uwm')) |>  
  left_join(trait_order) |>
  mutate(label = fct_reorder(label, -order), 
         group = fct_reorder(group, order), 
         term = as.character(term)) |> 
  filter(!is.na(order)) |> 
  group_by(variable) |> 
  filter(any(p.value < 0.05) & n > 4) |> 
  ungroup()

plot_change_pca <- plot_change_all4 |> 
  select(term, label, estimate) |> 
  pivot_wider(names_from = label, values_from = estimate) |> 
  mutate(across(everything(), ~ifelse(is.na(.x), 0.5, .x))) |> 
  as.data.frame() |> 
  column_to_rownames('term') 

pca <- rda(plot_change_pca~1, scale = T)
summary(plot_change_pca)

screeplot(pca, bstick = T)

eig <- paste0(names(eigenvals(pca)), ' (', round(eigenvals(pca)/sum(eigenvals(pca))*100, 1), '%)')

hab_scores <- scores(pca, choices = 1:3, display = 'sites', scaling = 'species') |> 
  as_tibble(rownames = 'habitat') |> 
  mutate(habitat = fct_relevel(habitat, 'A', 'M', 'R', 'TD_wet', 'TD_mesic', 'TE', 'TFG', 'TH', 'X'))

trait_scores <- scores(pca, choices = 1:3, display = 'species', scaling = 'species', const = 20) |> 
  as_tibble(rownames = 'trait') |> 
  mutate(fit_12 = goodness(pca, display = "sp", choices = 1:2, # goodness of fit at the 1st two axes
                           model = "CA", summarize = T)) 

trait_scores |> 
  ggplot(aes(PC1, PC2)) +
  # geom_vline(xintercept = 0, linetype = 'dashed', alpha = 0)+
  # geom_hline(yintercept = 0, linetype = 'dashed', alpha = 0)+
  geom_point(data = hab_scores, aes(fill = habitat), size = 5, pch = 21, alpha = 0.7)+
  geom_text_repel(aes(label = trait), max.overlaps = Inf, alpha = 0.7, color = '#6E2C00')+
  scale_fill_manual(values = veg_col, labels = c("Alpine and \nsubalpine vegetation",
                                                 "Wetland vegetation",
                                                 "Spring and mire \nvegetation",
                                                 "Wet meadows",
                                                 'Mesic meadows \nand pastures',
                                                 "Nardus grasslands \nand heathlands",
                                                 'Sand and shallow-soil \nvegetation',
                                                 'Dry grasslands',
                                                 "Ruderal and \nweed vegetation"))+
  guides(fill = guide_legend(byrow = T))+
  theme_bw()+
  theme(legend.position = 'bottom', legend.title = element_blank())+
  labs(x = eig[1], y = eig[2])

ggsave('plots/pca_traits_12.png', width = 8.5, height = 8.5)

