#' Supplementary code to the article: 
#' Klinkovska et al. Half a century of temperate non-forest vegetation changes: 
#' no net loss in species richness, but considerable shifts in taxonomic and functional composition.
#' 
#' Author: Klara Klinkovska, Helge Bruelheide 2024-09-24
#' R version 4.3.2

library(vegan) # version 2.6-4
library(broom) # version 1.0.5
library(mgcv) # version 1.9-1
library(tidyverse) # version 2.0.0
# dplyr version 1.1.3
# purrr version 1.0.2
# readr version 2.1.4
# tibble version 3.2.1
# tidyr version 1.3.1

# load data ---------------------------------------------------------------
traits <- read_csv('data/traits.csv')
spe <- read_csv('data/species.csv')
head <- read_csv('data/head.csv')

# plot characteristics ----------------------------------------------------

plot_traits <- spe |>
  left_join(traits) |> 
  mutate(Pres_abs = 1) |> 
  group_by(PlotObservationID) |> 
  summarise(richness = length(species_name), 
            shannon = diversity(cover_perc, 'shannon'), 
            evenness = shannon/log(richness), 
            across(EIV_light:log_height_mean, 
                   list(wm = function(x){weighted.mean(x, cover_perc, na.rm = T)}, 
                        uwm = function(x){weighted.mean(x, Pres_abs, na.rm = T)}), 
                   .names = "{.col}_{.fn}")) |> 
  write_csv('results/plot_characteristics.csv')

plot_traits <- read_csv('results/plot_characteristics.csv')

# extract trends - interval change approach -------------------------------------
head_int <- head |> 
  left_join(plot_traits) |> 
  pivot_longer(cols = richness:log_height_mean_uwm, names_to = 'variable', values_drop_na = T)

intervals <- head_int |> 
  arrange(desc(variable), RS_CODE_PLOT, year) |> 
  group_by(RS_CODE_PLOT, variable) |> 
  mutate(interval_change = ifelse(variable %in% c("richness", "shannon", "evenness"), 
                                  log(value/lag(value)), value - lag(value)), 
         from.PlotObservationID = lag(PlotObservationID), 
         to.PlotObservationID = PlotObservationID, 
         from = lag(year), to = year, 
         from_area = lag(Releve_area_m2), to_area = Releve_area_m2) |> 
  filter(!is.na(interval_change))

write_csv(intervals, 'results/intervals_all.csv')

# extract trends - linear change approach ---------------------------------------------------------

## ~ 5 min
trends <- plot_traits |> 
  pivot_longer(cols = richness:log_height_mean_uwm, names_to = 'variable') |> 
  filter(!is.na(value)) |> 
  left_join(head |>
              group_by(RS_CODE_PLOT) |> 
              mutate(min.year = min(year)) |> 
              ungroup() |> 
              mutate(year.fct = year - min.year) |> 
              select(PlotObservationID, RS_CODE_PLOT, year.fct, veg_type0), 
            by = "PlotObservationID") |> 
  nest(data = c(PlotObservationID, year.fct, value)) |> 
  mutate(m1 = map(data, ~lm(value~year.fct, data = .x)), 
         m_tidy = map(m1, ~tidy(.x)),
         est = map_dbl(m_tidy, ~.x$estimate[[2]]), 
         std.error = map_dbl(m_tidy, ~.x$std.error[[2]]),
         p = map_dbl(m_tidy, ~.x$p.value[[2]])) |> 
  select(RS_CODE_PLOT, veg_type0, variable, est, std.error, p)

write_csv(trends, 'results/lm_trends2.csv')
