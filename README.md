# Klinkovská et al. (2025) Global Change Biology: Half a century of temperate non-forest vegetation changes: no net loss in species richness, but considerable shifts in taxonomic and functional composition.

## Supplementary data and code to the article: 
Klinkovská K., Sperandii M. G., Knollová I., Danihelka J., Hájek M., Hájková P., Hroudová Z., Jiroušek M., Lepš J., Navrátilová J., Peterka T., Petřík P., Prach K., Řehounková K., Rohel J., Sobotka V., Vávra M., Bruelheide H. & Chytrý M. (2025) Half a century of temperate non-forest vegetation changes: no net loss in species richness, but considerable shifts in taxonomic and functional composition. Global Change Biology.

[![DOI](https://zenodo.org/badge/862873758.svg)](https://doi.org/10.5281/zenodo.14589270)

This repository contains vegetation-plot time series data and code used for the analysis of temporal trends in the Czech vegetation. 

## Data
The data contain plant species composition data from repeated vegetation-plot records from the Czech Republic. In total, the dataset consists of 1154 vegetation-plot time series from 53 resurvey studies comprising 3 909 vegetation-plot records sampled between 1971 and 2023. 

The data are provided as three CSV files with columns separated by commas:

* `Klinkovska_et_al_half_a_century_of_temperate_vegetation_change_species.csv` contains the percentage covers of plant species in the plots. Only records of vascular plants are included and the nomenclature and taxonomic concepts were standardized according to Kaplan et al. (2019). Most subspecies were merged to the species level and some species were merged into aggregates. For species records determined only to the genus level, we checked the source data, and if a species was determined at a lower taxonomic level in a different sampling event of the same plot, we related this record to the lower-level taxon (e.g., if Viola species was present in one time, and Viola hirta in another time in the same plot, Viola species was considered to be also Viola hirta). If more than one lower-level taxon occurred in another survey of the same plot, we equally distributed the cover of the genus-level record among the lower-level taxa. To minimize pseudoturnover caused by the misidentification of taxa in some surveys of a specific plot, we merged species we suspected to be misidentified under the name used in the last survey within a given time series. Moreover, we excluded vernal taxa Anemone nemorosa and Cardamine pratensis agg. from the plots in resurvey project CZ_0019_042 because the surveys were conducted in slightly different phenological stages (Klinkovská et al., 2023). We converted categories of different cover scales used to estimate species cover in vegetation plots to percentages representing the mean value of each interval. In some resurvey studies, different cover scales were used in the different surveys. In such cases, we converted the different cover scales into the least precise scale used in the time series (usually the nine-grade Braun-Blanquet scale to the seven-grade Braun-Blanquet scale; Westhoff & van der Maarel (1978)). 

* `Klinkovska_et_al_half_a_century_of_temperate_vegetation_change_head.csv` contains header data for the vegetation plots. The header data structure follows that of the ReSurveyEurope Database (http://euroveg.org/eva-database-re-survey-europe). 

* `Klinkovska_et_al_half_a_century_of_temperate_vegetation_change_traits.csv` contains species characteristics obtained from the Pladias Database of the Czech Flora and Vegetation (Chytrý et al., 2021). They include growth form (Dřevojan, 2020), life strategy scores (Guo and Pierce, (2019) following the method of Pierce et al., (2017)), height (Kaplan et al., 2019), leaf characteristics (E-Vojtkó et al., 2020; Findurová, 2018; Kleyer et al., 2008; Klotz & Kühn, 2002), flower characteristics (Durka, 2002), reproduction type (Chrtek, 2018; Durka, 2002), dispersal strategy (Sádlo et al., 2018), myrmecochory (Konečná et al., 2018), symbiosis with nitrogen fixers (Blažek & Lepš, 2016), trophic mode (Těšitel et al., 2016), taxon origin (Pyšek et al., 2022), Ellenberg-type indicator values (Chytrý et al., 2018), indicator values for disturbance of the herb layer (Herben et al., 2016), ecological specialization index (Zelený & Chytrý, 2019), indices of colonization ability (Prach et al., 2017), and Red List status (Grulich, 2017).

Original species composition data before nomenclature standardization and cover transformation and full header data are available in the ReSurveyEurope Database (http://euroveg.org/eva-database-re-survey-europe). 

## Scripts
* `1_extract_trends_community.R`: calculation of diversity indices, community weighted and unweighted means of functional and ecological species characteristics for each plot, calculation of trends in diversity indices and functional community characteristics using the interval change and linear trend approach. 
* `2_model_community.R`: generalized additive models for testing the probability of detecting a positive trend in each variable.
* `3_extract_trends_model_species.R`: changes at the species level.
* `4_species_composition.R`: distance-based redundancy analysis (db-RDA) for calculation of changes in species composition.
