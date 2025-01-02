#' Supplementary code to the article: 
#' Klinkovska et al. (2025) Half a century of temperate non-forest vegetation changes: 
#' no net loss in species richness, but considerable shifts in taxonomic and functional composition. Global Change Biology.
#' 
#' Author: Klara Klinkovska 2024-09-24
#' R version 4.3.2

library(vegan) # version 2.6-4
library(tidyverse) # version 2.0.0
# dplyr version 1.1.3
# readr version 2.1.4
# tibble version 3.2.1
# tidyr version 1.3.1

# load data ---------------------------------------------------------------
head <- read_csv('data/Klinkovska_et_al_half_a_century_of_temperate_vegetation_change_head.csv')

spe <- read_csv('data/Klinkovska_et_al_half_a_century_of_temperate_vegetation_change_species.csv') |> 
  pivot_wider(values_from = cover_perc, names_from = species_name, values_fill = 0)

### function for ordinations and results saving
rda_vegtype <- function(vegtype) {
  
  # filter one vegetation type
  print(vegtype)
  
  spe2 <- spe |> 
    semi_join(head |> 
                filter(veg_type0 == vegtype)) |> 
    arrange(PlotObservationID) |> 
    select(where(~sum(. > 0) > 1), 
           -c(PlotObservationID, RS_CODE)) |> 
    sqrt()
  
  head2 <- head |> 
    filter(veg_type0 == vegtype) |> 
    arrange(PlotObservationID)
  
  # constrained ordination - effect of time
  ord <- capscale(spe2 ~ head2$year + Condition(as.factor(head2$RS_CODE_PLOT)), 
                  distance = "bray", sqrt.dist = T)
  ord
  
  # significance test
  ord_anova <- anova(ord, permutations = how(blocks = as.factor(head2$RS_CODE_PLOT), nperm = 999)) |> 
    rownames_to_column('term')
  
  write_csv(ord_anova, paste0('results/ord_anova_', vegtype, '.csv'))
  
  eig <- paste0(names(eigenvals(ord)), ' (', round(eigenvals(ord)/ord$tot.chi*100, 1), '%)')[1:10] |> as_tibble()
  write_csv(eig, paste0('results/ord_eig_', vegtype, '.csv'))
}

for (i in unique(head$veg_type0)) {
  rda_vegtype(i)
}




