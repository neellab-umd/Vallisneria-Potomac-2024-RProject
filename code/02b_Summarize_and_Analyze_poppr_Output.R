## Author: Maile Neel
## 2024

## This script calculates all the genotypic diversity statistics in the manuscript except the distances among samples of the same MLG. It is designed to be run within the RProject called Potomac.Manuscript.2024.RProject.Rproj

#It uses the output of 02a_PR_Create_poppr_Files_and_Run_poppr.r and the original microsatellite data file. If that script has been run, you will find poppr.site.level.output.allpotomac.csv and poppr.tide.level.output.allpotomac.csv in the ./processed.data folder.  You can read those files in (code is below) to get the dataframes you need. If those files do not exist you need to run the script to create the necessary data objects.

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Load Packages & Source Functions ----
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,
               here,
               modelr,
               magrittr,
               Hmisc,
               broom,
               update = FALSE)

## Source custom functions needed for permutation testing.

if (!exists("summarize_by_tide.r", mode = "function"))
  source(here::here("functions", "summarize_by_tide.r"))

if (!exists("calculate_diffs", mode = "function"))
  source(here::here("functions", "permutation.functions.r"))

if (!exists("calculate_p_value", mode = "function"))
  source(here::here("functions", "permutation.functions.r"))


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Read data files into data.frames ----
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

overall.mlg.diversity <- readr::read_csv(here("processed.data",
                                       "poppr.tide.level.output.allpotomac.csv"))

mlg.diversity.by.site.with.labels<- readr::read_csv(here("processed.data",
                                                  "poppr.site.level.output.allpotomac.csv"))

gendata.allpotomac <- readr::read_csv(here("data",
                                    "allpotomac.microsatellite.data.csv")) %>% 
  drop_na(Clone.ID.2018)


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Summarize for the river and each tidal regime ----
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

## Get overall numbers of samples, MLG, and GD for the river and Tidal Regimes independent of sites from the overall poppr object created above


###  Number of MLGs ----

overall.mlg.diversity


### Number of Observations per MLG ----


# In this section we count the number of stems per MLG from the original source data file and calculate how many are single- versus multi-stemmed overall and in each tidal regime. We also get frequencies of each multisample clone by site for plotting in script 12_Mapping.With.ggplot.qmd

MLG.count <- gendata.allpotomac %>% 
  dplyr::count(Clone.ID.2018,Tide) %>% 
  dplyr::mutate(Category = ifelse(n > 1, "multi", "single")) %>% 
  dplyr::rename(stems = n)

## Number of MLGs with more than one observation versus with one observation by tidal group.
MultiSampleMLGs.by.Tide <- MLG.count %>% 
  bind_rows(., MLG.count %>% mutate(Tide = as.character("All"))) %>%
  dplyr::group_by(Tide, Category) %>%
  dplyr::summarize(Num.of.MLGs = n(), 
                   Num.of.Stems = sum(stems),
                   .groups = 'drop')


### Proportions of single vs. multistemmed individuals ----

#### Overall----

MLG.count %>% 
  dplyr::group_by(Category) %>%
  dplyr::summarize(Num.of.MLGs = n(), 
                   Num.of.Stems = sum(stems),
                   Proportion.of.MLGs = n()/length(unique(MLG.count$Clone.ID.2018)),
                   Proportion.of.stems = sum(stems)/sum(MLG.count$stems),
                   .groups = 'drop')

#### Within the Tidal----

number.tidal.MLGs <-MLG.count %>% 
  dplyr::filter(Tide == "Tidal") %>% dplyr::summarize(length(unique(Clone.ID.2018))) %>% 
  as.vector()

number.tidal.stems <- MLG.count %>% 
  dplyr::filter(Tide == "Tidal") %>% dplyr::summarize(sum(stems)) %>% 
  as.vector()

MLG.count %>% 
  dplyr::filter(Tide == "Tidal") %>% 
  dplyr::group_by(Category) %>%
  dplyr::summarize(Num.of.MLGs = n(), 
                   Num.of.Stems = sum(stems),
                   Proportion.of.MLGs = n()/number.tidal.MLGs[[1]],
                   Proportion.of.stems = sum(stems)/number.tidal.stems[[1]],
                   .groups = 'drop')

#### Within the NonTidal -----

number.nontidal.MLGs <-MLG.count %>% 
  dplyr::filter(Tide == "NonTidal") %>% dplyr::summarize(length(unique(Clone.ID.2018))) %>% 
  as.vector()

number.nontidal.stems <- MLG.count %>% 
  dplyr::filter(Tide == "NonTidal") %>% dplyr::summarize(sum(stems)) %>% 
  as.vector()

MLG.count %>% 
  dplyr::filter(Tide == "NonTidal") %>% 
  dplyr::group_by(Category) %>%
  dplyr::summarize(Num.of.MLGs = n(), 
                   Num.of.Stems = sum(stems),
                   Proportion.of.MLGs = n()/number.nontidal.MLGs[[1]],
                   Proportion.of.stems = sum(stems)/number.nontidal.stems[[1]],
                   .groups = 'drop')

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Summarize Multisample MLGS by site ----
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

## Total number of stems sampled at each site ----

# The total number of stems is needed to calculated frequencies of each multisample MLG within each site.
Stem.Count.By.Site <-  gendata.allpotomac %>% 
  dplyr::group_by(NewPop,OrderPop, Tide) %>% 
  dplyr::count() %>% 
  dplyr::rename(Number_of_Stems = n)

## Counts and Frequencies of MLGs ----

## Get site-level counts and calculate frequencies of all multisample MLGs sampled at multiple nontidal sites
MultiSampleMLG.Count.by.Site <- gendata.allpotomac %>% 
  dplyr::group_by(Clone.ID.2018) %>% # get all instances of MLGS that occur more than once
  dplyr::filter(n()>1) %>% 
  dplyr::ungroup() %>% ## get other variables for each sample back into the data.frame
  dplyr::group_by(NewPop) %>% 
  dplyr::count(Clone.ID.2018) %>% 
  dplyr::rename(MLG = n) %>% 
  dplyr::group_by(Clone.ID.2018) %>% # get all instances of MLGS that occur at >1 site
  dplyr::filter(n()>1) %>% 
  dplyr::ungroup() %>% ## get other variables for each sample back into the data.frame
  full_join(Stem.Count.By.Site, by = "NewPop") %>%
  filter(Tide == "NonTidal") %>% 
  dplyr::mutate(FR = MLG/Number_of_Stems) %>% 
  dplyr::arrange(Clone.ID.2018) %>% 
  tidyr::pivot_wider(names_from = c(Clone.ID.2018),values_from = c(MLG, FR), values_fill = 0)

## Write out file for later plotting ----
write_csv(MultiSampleMLG.Count.by.Site, here("processed.data", "mlg.count.and.freq.by.pop.csv" ))


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Summary statistics for MLGs Within Sites Overall and By Tidal Regime ----
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

### Site level genotypic diversity overall and by tidal regime ----

mlg.diversity.by.site.with.labels %>% 
  dplyr::select(NewPop, Tide, N:MLG,Eff., SW:GD.Simp) %>% 
  summarize_by_tide()


### Rarefied estimates of site level genotypic diversity overall and by tidal regime ----

mlg.diversity.by.site.with.labels %>% 
  dplyr::select(NewPop, Tide, eMLG, Eff.SW.est.5:GD.Simp.est.5) %>% 
  summarize_by_tide()


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Test for differences using permutation tests ----
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

## permute creates the chosen number of dataframes in each of which each chosen variable has been permuted. This action creates a list in which each element is a permuted data frame.

permutations <- modelr::permute(mlg.diversity.by.site.with.labels, 10000, MLG:GD.Simp.est.5)

## Use our calculate_diffs() function to calculate differences between means of values in NonTidal versus Tidal environments.   The variable Tide is hard-wired into the calculate_diffs() function. If you want to permute by a different variable you need to modify the function.

### Do the calculation on the  observed values
observed_diffs<-mlg.diversity.by.site.with.labels %>% 
  dplyr::select(-c(OrderPop, N, Eff.SW.Lower.CI.5:Eff.Simp.Upper.CI.5)) %>% 
calculate_diffs()

## map our function over all elements in the list that was output by permute() to test the differences between the means in Tidal versus NonTidal environments on all the permuted values - this step can be time consuming if you have lots of permutations

permuted_diffs_list <- purrr::map(permutations$perm, calculate_diffs) 

## Pull differences from permuted data out of the list and into a data_frame
permuted_diffs_df <-do.call("rbind", permuted_diffs_list) 

# map our calculate_p_value() function over the distribution of permuted differences to calculate the number of permutations that have a more extreme difference than our observed value. this number out of the number of permutations is the p-value.

p_values_tibble <- map2(names(observed_diffs), 
                        observed_diffs, 
                        ~ {calculate_p_value(.x, .y, permuted_diffs_df[[.x]])
                        }) %>%
  bind_rows()

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#  Correlations with River Km ----
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

## Read in river km and join that variable with mlg data.

RiverKm <- readr::read_csv(here("data","PR.Pop.Centroids.UTMS.csv")) %>% 
  dplyr:: select(NewPop, RiverKm)

mlg.diversity.by.site.with.labels <- mlg.diversity.by.site.with.labels %>% 
  dplyr::left_join(.,RiverKm, by= "NewPop")

## Function for running correlations on multiple variables

corrfunct <- function(df,vars) {
  df %>% 
    dplyr::select(all_of(vars), RiverKm) %>% 
    as.matrix() %>% 
    Hmisc::rcorr(type = "spearman")
}

## Choose variables for correlations ----

vars.mlg <-c("MLG", "eMLG", "Eff.SW", "Eff.Simp","GD", "GD.SW","GD.Simp", "Eff.SW.est.5","Eff.Simp.est.5","GD.SW.est.5", "GD.Simp.est.5")

## Correlations - Full River ----

corr.results.fullriver.mlg <- map(vars.mlg, 
                            corrfunct, 
                            df = mlg.diversity.by.site.with.labels) %>% 
  purrr::set_names(.,vars.mlg) %>% 
  purrr::map_dfr(broom::tidy) %>% 
  dplyr::select(-column1, -n) %>% 
  dplyr::rename(r.full.river = estimate, p.full.river = p.value)

## Correlations - Nontidal ----

corr.results.Tidal.mlg <- map(vars.mlg, 
                          corrfunct, 
                          df=filter(mlg.diversity.by.site.with.labels, 
                                    Tide=="Tidal")) %>% 
  purrr::set_names(.,vars.mlg) %>% 
  purrr::map_dfr(broom::tidy) %>% 
  dplyr::select(-column1, -n) %>% 
  rename(r.tidal = estimate, p.tidal = p.value)

## Correlations - Tidal ----

corr.results.Non.Tidal.mlg <- purrr::map(vars.mlg, 
                              corrfunct, 
                              df=filter(mlg.diversity.by.site.with.labels, 
                                        Tide == "NonTidal")) %>% 
  purrr::set_names(.,vars.mlg) %>% 
  purrr::map_dfr(broom::tidy) %>% 
  dplyr::select(-column1, -n) %>% 
  dplyr::rename(r.nontidal = estimate, p.nontidal = p.value)

## Correlations - Summary ----

correlation.table.mlg<-corr.results.fullriver.mlg %>% 
  dplyr::left_join(corr.results.Non.Tidal.mlg, by = "column2") %>% 
  dplyr::left_join(corr.results.Tidal.mlg, by = "column2")

