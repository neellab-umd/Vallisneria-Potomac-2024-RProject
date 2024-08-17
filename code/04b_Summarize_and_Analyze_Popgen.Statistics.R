## Author: Maile Neel
## 2024

#This script os designed to be used within the Potomac.Manuscript.2024.R.Project.Rproj RProject. It calculates summaries for raw population genetic statistics stored in objects that were created in the script 04a_Calculate.Popgen.Statistics.R.  The script summarize to the overall river and by tidal regime.  It runs permutation tests to evaulate seeing diversity levels in tidal sites versus nontidal sites by chance. 

## These required objects were written to csv files in the ./processed data folder.

## final.all.global.among.population.measures.csv
## FST_confidence.intervals.csv 
## MMOD_confidence.intervals.df.csv" 
## final.site.level.popgen.stats.with.labels.csv 

## If the objects containing the population genetic statistics are not in the environment or saved in the processed.data folder, that script needs to be run first. 

## RiverKm from the PR.Pop.Centroids.UTMs.csv from the ./data folder is also used.

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Install and Load Packages & Source Functions ----
#XXXXXXXXXXXXXXXXXXXXxxxXXXXXXXXXXXXXXXXXXXXXXXXXXXX

if (!require("pacman")) install.packages("pacman")
pacman::p_load(here,
               tidyverse,
               modelr,
               Hmisc,
               broom,
               boot)

## Source custom functions needed for permutation testing.

if (!exists("calculate_diffs", mode = "function"))
  source(here::here("functions", "permutation.functions.r"))

if (!exists("calculate_p_value", mode = "function"))
  source(here::here("functions", "permutation.functions.r"))

if (!exists("sumarize_by_tide", mode = "function"))
  source(here::here("functions", "summarize_by_tide.r"))


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
## Read in csv files  of global among population diversity estimates ----
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

final.all.global.among.population.measures <- readr::read_csv(here("processed.data", 
                                                                   "final.all.global.among.population.measures.csv"))

#Get Confidence intervals from the objects created or read them back in if they are not in the environment.

MMOD_confidence.intervals.df <- readr::read_csv(here("processed.data",
                                                     "MMOD_confidence.intervals.df.csv"))

FST_confidence.intervals <- readr::read_csv(here("processed.data",
                                                 "FST_confidence.intervals.csv"))

#+++++++++++++++++++++++++++++++++++++++++++++++
## Summarize Pairwise Values Among Sites ----
#+++++++++++++++++++++++++++++++++++++++++++++++

#read in file with summary values for site-level statistics if the object is not already in the environment.

final.site.level.popgen.stats.with.labels <- 
  readr::read_csv( here("processed.data", 
                        "final.site.level.popgen.stats.with.labels.csv"))

## For ease of getting numbers from the output for the manuscript text, the summaries are done in two smaller groups

## Within Site Statistics

final.site.level.popgen.stats.with.labels %>% 
  dplyr::select(c(NewPop,Tide,PercPolymorphic:Fis)) %>% 
  summarize_by_tide()

##Among-Site, River-wide & within-Tidal Regime

final.site.level.popgen.stats.with.labels %>% 
  dplyr::select(c(NewPop,Tide,FST:Josts.D.within)) %>% 
  summarize_by_tide()

#+++++++++++++++++++++++++++++++++++++++++++
## Permutation tests for differences ----
#+++++++++++++++++++++++++++++++++++++++++++

## permute creates the chosen number of dataframes in each of which each chosen variable has been permuted. This action creates a list in which each element is a permuted data frame.

permutations.pop.gen <- modelr::permute(final.site.level.popgen.stats.with.labels, 10000, c(PercPolymorphic:Fis,FST:Josts.D.within))

## Use our calculate_diffs() function to calculate differences between means of values in NonTidal versus Tidal environments.   The variable Tide is hard-wired into the calculate_diffs() function. If you want to permute by a different variable you need to modify the function.

### Do the calculation on the  observed values
observed_diffs.pop.gen <- final.site.level.popgen.stats.with.labels %>% 
  dplyr::select(-c(OrderPop, N.Potomac.No.Reps, Lower_Fis_CI, Upper_Fis_CI)) %>% 
  calculate_diffs()

## map our function over all elements in the list that was output by permute() to test the differences between the means in Tidal versus NonTidal environments on all the permuted values - this step can be time consuming if you have lots of permutations

permuted_diffs_list.pop.gen <- map(permutations.pop.gen$perm, calculate_diffs) 

## Pull differences from permuted data out of the list and into a data_frame
permuted_diffs_df.pop.gen <-do.call("rbind", permuted_diffs_list.pop.gen) 

# map our calculate_p_value() function over the distribution of permuted differences to calculate the number of permutations that have a more extreme difference than our observed value. this number out of the number of permutations is the p-value.

p_values.pop.gen.tibble <- purrr::map2(names(observed_diffs.pop.gen), 
                                observed_diffs.pop.gen, 
                        ~ {calculate_p_value(.x, .y, permuted_diffs_df.pop.gen[[.x]])
                        }) %>%
  bind_rows()

#++++++++++++++++++++++++++++++++++
##  Correlations with River Km ----
#++++++++++++++++++++++++++++++++++

## Read in river km and join that variable with mlg data.

RiverKm <- readr::read_csv(here("data","PR.Pop.Centroids.UTMS.csv")) %>% 
  dplyr:: select(NewPop, RiverKm)

final.site.level.popgen.stats.with.labels <- final.site.level.popgen.stats.with.labels %>% 
left_join(.,RiverKm, by= "NewPop")

## Function for running correlations on multiple variables

corrfunct <- function(df,vars) {
df %>% 
  dplyr::select(all_of(vars),RiverKm) %>% 
  as.matrix() %>% 
  Hmisc::rcorr(type="spearman")
}

## Choose variables for correlations ----

vars.pgen <-c("Ap","Ar","Ho","He","Fis")

## Correlations - Full River ----

corr.results.fullriver<-map(vars.pgen, 
                            corrfunct, 
                            df=final.site.level.popgen.stats.with.labels) %>% 
  purrr::set_names(.,vars.pgen) %>% 
  purrr::map_dfr(broom::tidy) %>% 
  dplyr::select(-column1, -n) %>% 
  dplyr::rename(r.full.river=estimate, p.full.river =p.value)

## Correlations - Tidal ----

corr.results.Tidal <- map(vars.pgen, 
                            corrfunct, 
                            df=filter(final.site.level.popgen.stats.with.labels, 
                                      Tide=="Tidal")) %>% 
  purrr::set_names(.,vars.pgen) %>% 
  purrr:: map_dfr(broom::tidy) %>% 
  dplyr::select(-column1, -n) %>% 
  dplyr::rename(r.tidal=estimate, p.tidal =p.value)

## Correlations - Nontidal ----

corr.results.Non.Tidal <- map(vars.pgen, 
                          corrfunct, 
                          df=filter(final.site.level.popgen.stats.with.labels, 
                                    Tide=="NonTidal")) %>% 
  purrr::set_names(.,vars.pgen) %>% 
  purrr::map_dfr(broom::tidy) %>% 
  dplyr::select(-column1, -n) %>% 
  dplyr::rename(r.nontidal=estimate, p.nontidal =p.value)

## Correlations - Summary ----

correlation.table<-corr.results.fullriver %>% 
  left_join(corr.results.Non.Tidal, by="column2") %>% 
  left_join(corr.results.Tidal, by="column2")

