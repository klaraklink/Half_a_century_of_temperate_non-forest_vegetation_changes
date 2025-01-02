#' Supplementary code to the article: 
#' Klinkovska et al. (2025) Half a century of temperate non-forest vegetation changes: 
#' no net loss in species richness, but considerable shifts in taxonomic and functional composition. Global Change Biology.
#' 
#' Author: Klara Klinkovska, Helge Bruelheide 2024-09-24
#' R version 4.3.2

library(broom) # version 1.0.5
library(sf) # version 1.0-14
library(patchwork) # version 1.2.0
library(tidyverse) # version 2.0.0
# dplyr version 1.1.3
# ggplot2 version 3.5.1
# purrr version 1.0.2
# readr version 2.1.4
# stringr version 1.5.1
# tibble version 3.2.1
# tidyr version 1.3.1


# load data ---------------------------------------------------------------

head <- read_csv('data/Klinkovska_et_al_half_a_century_of_temperate_vegetation_change_head.csv')

spe <- read_csv('data/Klinkovska_et_al_half_a_century_of_temperate_vegetation_change_species.csv') |> 
  left_join(head |> 
              select(PlotObservationID, RS_CODE_PLOT, year, Releve_area_m2, veg_type0) |> 
              group_by(RS_CODE_PLOT) |> 
              mutate(min.year = min(year)) |> 
              ungroup() |> 
              mutate(year.fct = year - min.year, .after = min.year)) |> 
  anti_join(read_csv('data/diff_area.csv') |> 
  select(RS_CODE_PLOT = RS_CODE_PLOT.x, from.PlotObservationID, to.PlotObservationID, 
         Releve_area_m2.x, Releve_area_m2.y))

# join shp for grid to select only species with records from more grid cells
head_grid <- head |> 
  group_by(RS_CODE_PLOT) |> 
  slice_max(year) |> 
  ungroup() |> 
  st_as_sf(coords = c('Longitude', 'Latitude'), crs = 4326) |> 
  st_join(read_sf('shps/pladias_qq4_CZ_2017_WGS84.shp') |> 
            select(code, square)) |> 
  arrange(PlotObservationID) |> 
  as_tibble()


# extract trends (lm) for species  --------------------------------------------------

# absolute changes in cover
spe_lm <- spe |> 
  pivot_wider(names_from = species_name, values_from = cover_perc, values_fill = 0) |> 
  pivot_longer(cols = -c(1:8), names_to = 'species_name', values_to = 'cover_perc') |> 
  arrange(RS_CODE_PLOT, species_name, year) |> 
  group_by(RS_CODE_PLOT, species_name) |> 
  mutate(cover_all = sum(cover_perc), year) |> 
  filter(cover_all != 0) |> 
  nest(data = -c(RS_CODE_PLOT, species_name, veg_type0)) |> 
  mutate(m1 = map(data, ~lm(cover_perc~year.fct, data = .x)), 
         m_tidy = map(m1, ~tidy(.x)), 
         lm_change = map_dbl(m_tidy, ~.x$estimate[[2]]), 
         p = map_dbl(m_tidy, ~.x$p.value[[2]])) |> 
  select(RS_CODE_PLOT, veg_type0, species_name, lm_change, p)

write_csv(spe_lm, 'results/lm_species.csv')
spe_lm <- read_csv('results/lm_species.csv')

# absolute changes p/a 
spe_lm_pa <- spe |> 
  #filter(species_name %in% c('Agrimonia eupatoria', 'Arrhenatherum elatius')) |> #test
  mutate(Pres = ifelse(cover_perc > 0, 1, 0), .before = 'cover_perc', .keep = 'unused') |> # replace cover values by presence absence data
  pivot_wider(names_from = species_name, values_from = Pres, values_fill = 0) |> # add zeros
  pivot_longer(cols = -c(1:8), names_to = 'species_name', values_to = 'Pres') |> 
  arrange(RS_CODE_PLOT, species_name, year) |> 
  group_by(RS_CODE_PLOT, species_name) |> 
  mutate(pres_all = sum(Pres)) |> 
  filter(!all(Pres == 0) & !all(Pres == 1)) |> 
  nest(data = -c(RS_CODE_PLOT, species_name, veg_type0)) |> 
  mutate(m1 = map(data, ~glm(Pres~year.fct, family = 'binomial', data = .x)), 
         m_tidy = map(m1, ~tidy(.x)), 
         lm_change = map_dbl(m_tidy, ~.x$estimate[[2]]), 
         p = map_dbl(m_tidy, ~.x$p.value[[2]])) |> 
  select(RS_CODE_PLOT, veg_type0, species_name, lm_change, p)

write_csv(spe_lm_pa, 'results/lm_species_pa.csv')

# join together cover and pa 
spe_lm_all <- bind_rows(read_csv('results/lm_species.csv') |> 
                          mutate(change = 'cover'), 
                        read_csv('results/lm_species_pa.csv') |> 
                          mutate(change = 'pa')) |> 
  filter(str_detect(species_name, 'species', negate = T) | 
           species_name %in% c('Crataegus species', 'Oenothera species'))

# check
spe_lm_all  |> 
  distinct(species_name)

# linear trends proportion of increasing, decreasing trends -------------------------------------------
# binomial test
spe_lm2 <- spe_lm_all |> 
  filter(!(change == 'pa' & near(lm_change, 0))) |> 
  left_join(head_grid) |> 
  group_by(species_name, change) |> 
  summarise(n = n(),
            n_grid = n_distinct(code),
            n_pos = length(lm_change[lm_change > 0]), 
            n_neg = length(lm_change[lm_change < 0]), 
            n_all = ifelse(n_pos + n_neg > 0, n_pos + n_neg, 1), 
            est = binom.test(n_pos, n_all)$est, 
            p.values.binom = binom.test(n_pos, n_all)$p.value, 
            conf.binom.minus =  binom.test(n_pos, n_all)$conf.int[1],
            conf.binom.plus =  binom.test(n_pos, n_all)$conf.int[2],
            mean.absolute.change = mean(lm_change, na.rm = T),
            .groups = 'drop') 

# number of significantly increasing, decreasing species
spe_lm2 |> 
  filter(p.values.binom < 0.05) |> 
  group_by(sign(-(0.5 - est)), change) |> 
  count()

spe_lm2 |> 
  group_by(sign(mean.absolute.change), change) |> 
  count()


# plots -------------------------------------------------------------------

# all significantly increasing and decreasing species
spe_lm2 |> 
  filter(change == 'pa' & p.values.binom < 0.05) |>
  arrange(est) |>
  mutate(ID = row_number(), species_name = factor(species_name, ordered = T, level = species_name)) |> 
  ggplot(aes(-(0.5 - est), species_name)) +
  geom_col(fill = 'green4', color = 1, alpha = 0.5)+
  geom_errorbar(aes(xmin = -(0.5 - conf.binom.minus), xmax = -(0.5 - conf.binom.plus)), width = 0.5)+
  scale_x_continuous(breaks = seq(-0.5, 0.5, 0.25), labels = seq(0, 1, 0.25))+
  facet_wrap(~change, labeller = labeller(change = c(pa = '(a) Presence-absence')))+
  theme_bw()+
  theme(axis.title.y = element_blank(), strip.background = element_blank())+
  labs(x = 'Probability to increase')+
  
  spe_lm2 |>  
  filter(change == 'cover' & p.values.binom < 0.05) |> 
  arrange(est) |> 
  mutate(ID = row_number(), species_name = factor(species_name, ordered = T, level = species_name)) |> 
  ggplot(aes(-(0.5 - est), species_name)) +
  geom_col(fill = 'green4', color = 1, alpha = 0.5)+
  geom_errorbar(aes(xmin = -(0.5 - conf.binom.minus), xmax = -(0.5 - conf.binom.plus)), width = 0.5)+
  scale_x_continuous(breaks = seq(-0.5, 0.5, 0.25), labels = seq(0, 1, 0.25))+
  facet_wrap(~change, labeller = labeller(change = c(cover = '(b) Cover')))+
  theme_bw()+
  theme(axis.title.y = element_blank(), strip.background = element_blank())+
  labs(x = 'Probability to increase in cover')

ggsave('plots/species_binom_lm_cover_pa.png', width = 9, height = 25)

# subset for the main article
# min number of observations of change + min number of sites = grid cells
spe_lm2 |> 
  filter(change == 'pa' & p.values.binom < 0.05 & n_all >= 30 & n_grid >= 15) |>
  arrange(est) |>
  mutate(ID = row_number(), species_name = factor(species_name, ordered = T, level = species_name)) |> 
  ggplot(aes(-(0.5 - est), species_name)) +
  geom_col(fill = 'green4', color = 1, alpha = 0.5)+
  geom_errorbar(aes(xmin = -(0.5 - conf.binom.minus), xmax = -(0.5 - conf.binom.plus)), width = 0.5)+
  scale_x_continuous(breaks = seq(-0.5, 0.5, 0.25), labels = seq(0, 1, 0.25))+
  facet_wrap(~change, labeller = labeller(change = c(pa = '(a) Presence-absence')))+
  theme_bw()+
  theme(axis.title.y = element_blank(), strip.background = element_blank())+
  labs(x = 'Probability to increase')+
  
  spe_lm2 |>  
  filter(change == 'cover' & p.values.binom < 0.05 & n_all >= 30 & n_grid >= 15) |> 
  arrange(est) |> 
  mutate(ID = row_number(), species_name = factor(species_name, ordered = T, level = species_name)) |> 
  ggplot(aes(-(0.5 - est), species_name)) +
  geom_col(fill = 'green4', color = 1, alpha = 0.5)+
  geom_errorbar(aes(xmin = -(0.5 - conf.binom.minus), xmax = -(0.5 - conf.binom.plus)), width = 0.5)+
  scale_x_continuous(breaks = seq(-0.5, 0.5, 0.25), labels = seq(0, 1, 0.25))+
  facet_wrap(~change, labeller = labeller(change = c(cover = '(b) Cover')))+
  theme_bw()+
  theme(axis.title.y = element_blank(), strip.background = element_blank())+
  labs(x = 'Probability to increase in cover')

ggsave('plots/species_binom_lm_cover_pa_subset.png', width = 8, height = 11.5)


# species by habitats -----------------------------------------------------
# color settings
veg_col <- c('A' = "#F8766D", 'M' = "#D39200", 'R' = "#DB72FB", 'TD_mesic' = "#00B9E3", 'TD_wet' = "#619CFF",
             'TE' = "#00C19F", 'TFG' = "#00BA38", 'TH' = "#A3A500", 'X' = "#FF61CC")

spe_lm3 <- spe_lm_all |> 
  group_by(species_name, veg_type0, change) |> 
  summarise(n = n(),
            n_pos = length(lm_change[lm_change > 0]), 
            n_neg = length(lm_change[lm_change < 0]), 
            n_all = ifelse(n_pos + n_neg > 0, n_pos + n_neg, 1), 
            est = binom.test(n_pos, n_all)$est, 
            p.values.binom = binom.test(n_pos, n_all)$p.value, 
            conf.binom.minus =  binom.test(n_pos, n_all)$conf.int[1],
            conf.binom.plus =  binom.test(n_pos, n_all)$conf.int[2],
            mean.absolute.change = mean(lm_change, na.rm = T), .groups = 'drop')

# numbers of significantly increasing, decreasing species
incr_decr <- spe_lm3 |> 
  group_by(sign(-(0.5 - est)), change, veg_type0) |> 
  count()

# ggplot for the most increasing and decreasing species
change_plot <- spe_lm3 |>  
  arrange(est, species_name) |> 
  filter(p.values.binom < 0.05 & sign(-(0.5 - est)) == sign(mean.absolute.change)) |> 
  mutate(rowID = row_number()) 

# plot species in one habitat type
spe_plot_hab <- function(vegtype, cover_pa) {
  change_plot |> 
    filter(veg_type0 == vegtype & change == cover_pa) |> 
    mutate(species_name = factor(species_name, labels = species_name, levels = species_name)) |> 
    ggplot(aes(-(0.5 - est), species_name, fill = veg_type0)) +
    geom_col(color = 1, alpha = 0.5)+
    geom_errorbar(aes(xmin = -(0.5 - conf.binom.minus), xmax = -(0.5 - conf.binom.plus)), width = 0.5)+
    scale_x_continuous(breaks = seq(-0.5, 0.5, 0.25), labels = seq(0, 1, 0.25))+
    scale_fill_manual(values = veg_col)+
    facet_wrap(~change, scales = 'free', 
               labeller = labeller(change = c(pa = '(a) Presence-absence', 
                                              cover = '(b) Cover'))) +
    theme_bw()+
    theme(axis.title.y = element_blank(), legend.position = 'none', strip.background = element_blank(), 
          strip.text = element_text(size = 12))
}

spe_plot_hab('A', 'pa') +
  labs(x = 'Probability to increase') + 
  spe_plot_hab('A', 'cover') +
  labs(x = 'Probability to increase in cover') +
  plot_annotation('Alpine and subalpine vegetation') 

ggsave('plots/species_binom_lm_cover_pa_A.png', width = 8, height = 3)

spe_plot_hab('M', 'pa') +
  labs(x = 'Probability to increase') + 
  spe_plot_hab('M', 'cover') +
  labs(x = 'Probability to increase in cover') +
  plot_annotation('Wetland vegetation')

ggsave('plots/species_binom_lm_cover_pa_M.png', width = 8, height = 4)

spe_plot_hab('R', 'pa') +
  labs(x = 'Probability to increase') + 
  spe_plot_hab('R', 'cover') +
  labs(x = 'Probability to increase in cover') +
  plot_annotation('Spring and mire vegetation')

ggsave('plots/species_binom_lm_cover_pa_R.png', width = 8, height = 4)

spe_plot_hab('TD_wet', 'pa') +
  labs(x = 'Probability to increase') + 
  spe_plot_hab('TD_wet', 'cover') +
  labs(x = 'Probability to increase in cover') +
  plot_annotation('Wet meadows')

ggsave('plots/species_binom_lm_cover_pa_TD_wet.png', width = 8, height = 8)

spe_plot_hab('TD_mesic', 'pa') +
  labs(x = 'Probability to increase') + 
  spe_plot_hab('TD_mesic', 'cover') +
  labs(x = 'Probability to increase in cover') +
  plot_annotation('Mesic meadows and pastures')

ggsave('plots/species_binom_lm_cover_pa_TD_mesic.png', width = 8, height = 3)

spe_plot_hab('TE', 'pa') +
  labs(x = 'Probability to increase') + 
  spe_plot_hab('TE', 'cover') +
  labs(x = 'Probability to increase in cover') +
  plot_annotation('Nardus grasslands and heathlands')

ggsave('plots/species_binom_lm_cover_pa_TE.png', width = 8, height = 3)

spe_plot_hab('TFG', 'pa') +
  labs(x = 'Probability to increase') + 
  spe_plot_hab('TFG', 'cover') +
  labs(x = 'Probability to increase in cover') +
  plot_annotation('Sand and shallow-soil vegetation')

ggsave('plots/species_binom_lm_cover_pa_TFG.png', width = 8, height = 4)

spe_plot_hab('TH', 'pa') +
  labs(x = 'Probability to increase') + 
  spe_plot_hab('TH', 'cover') +
  labs(x = 'Probability to increase in cover') +
  plot_annotation('Dry grasslands')

ggsave('plots/species_binom_lm_cover_pa_TH.png', width = 9, height = 14)

spe_plot_hab('X', 'pa') +
  labs(x = 'Probability to increase') + 
  spe_plot_hab('X', 'cover') +
  labs(x = 'Probability to increase in cover') +
  plot_annotation('Ruderal and weed vegetation')

ggsave('plots/species_binom_lm_cover_pa_X.png', width = 8, height = 6)
