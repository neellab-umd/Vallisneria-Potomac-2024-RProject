## Author: Maile Neel
## 2024

# Usage ----

## This script does the permutation tests to determine if relatedness differs in the different tidal environements and whether relatedness is correlated with river kilometers.  It is designed to work in the Potomac.Manuscript.2024.R.Project.RProj.  

## It runs permutation and correlation tests on summary dataframes of Wang's r created using the related package. Wang's r was calculated in script 06a_Calculate_Relatedness.R.  The summaries were done in 06b_Summarize_Relatedness_Output.R. The files saved in that script were written to the related subfolder in the processed.data folder in this project. 

## Files for all samples are in ./processed.data/related/allsamples, and files for distinct MLGs are in ./processed.data/related/noreps. If the files are not present in those folders, the 06a_Calculate_Relatedness.R and/or. 06b_Summarize_Relatedness_Output.R script need to be run to create them.

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# INSTALL AND LOAD PACKAGES ####
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,
               here,
               Hmisc,
               modelr
               )

#XXXXXXXXXXXXXXXXXXXXXX
# SOURCE FUNCTIONS ####
#XXXXXXXXXXXXXXXXXXXXXX

## Source custom functions needed for permutation testing.

if (!exists("calculate_diffs", mode = "function"))
  source(here::here("functions", "permutation.functions.r"))

if (!exists("calculate_p_value", mode = "function"))
  source(here::here("functions", "permutation.functions.r"))

## Source custom function needed to calculated means by tide and overall.

if (!exists("sumarize_by_tide", mode = "function"))
  source(here::here("functions", "summarize_by_tide.r"))

## Source function for running correlations with river km

corrfunct <- function(df,vars) {
  df %>% 
    dplyr::select(all_of(vars),RiverKm) %>% 
    as.matrix() %>% 
    Hmisc::rcorr(type="spearman")
}


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#         All SAMPLES      ####
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#+++++++++++++++++++
### Read in data ####
#+++++++++++++++++++

## Read in all Wangs R for all samples. Includes rw and rA for sample-level comparisons. Take out all sample pairs that are from different tidal regimes.

raw_relatedness_allsamples_wang_bothways.sametide <- 
  readr::read_csv(here("processed.data",
                       "related","allsamples", 
        "raw_relatedness_allsamples_wang_bothways.csv")) %>% 
  filter(Tide.Comp == "Same.Tide") 


## Read in rW and wA for all samples.
rW.and.rA.allsamples.by.site <- readr::read_csv("./processed.data/related/allsamples/wang.rW.and.rA.allsamples.by.site.csv")

#++++++++++++++++++++
## Permutations ####
#++++++++++++++++++++

#++++++++++++++++++++++++++++++++++++++++++++++++++++
#### Permute Overall r between Tidal and NonTidal  ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++

## permute creates the chosen number of dataframes in each of which each chosen variable has been permuted. This action creates a list in which each element is a permuted data frame.

permutations.all <- modelr::permute(raw_relatedness_allsamples_wang_bothways.sametide, 10000, Wang)

## Use our calculate_diffs() function to calculate differences between means of values in NonTidal versus Tidal environments.  The variable Tide is hard-wired into the calculate_diffs() function. If you want to permute by a different variable you need to modify the function.

### Do the calculation on the  observed values
observed_diffs.wang.all<-raw_relatedness_allsamples_wang_bothways.sametide %>% 
  dplyr::select(c(Tide,Wang)) %>% 
  calculate_diffs()

## map our function over all elements in the list that was output by permute() to test the differences between the means in Tidal versus NonTidal environments on all the permuted values - this step can be time consuming if you have lots of permutations

permuted_diffs_list.wang.all <- map(permutations.all$perm, calculate_diffs) 

## Pull differences from permuted data out of the list and into a data_frame
permuted_diffs_df.wang.all <-do.call("rbind", permuted_diffs_list.wang.all) 

# map our calculate_p_value() function over the distribution of permuted differences to calculate the number of permutations that have a more extreme difference than our observed value. this number out of the number of permutations is the p-value.

p_values_tibble.wang.all <- purrr::map2(names(observed_diffs.wang.all), 
                        observed_diffs.wang.all, 
                        ~ {calculate_p_value(.x, .y, permuted_diffs_df.wang.all[[.x]])
                        }) %>%
  bind_rows()

#+++++++++++++++++++++++++++++++++++
#### Permute rA for all samples.  #### 
#+++++++++++++++++++++++++++++++++++

## Set up the permutation.

permutations.wangtype.rA.allsamples <- 
  raw_relatedness_allsamples_wang_bothways.sametide %>% 
  dplyr::filter(wangtype == "rA") %>% 
  modelr::permute(., 10000, Wang)

## Use our calculate_diffs() function to calculate observed differences between means of values in NonTidal versus Tidal environments.  The variable Tide is hard-wired into the calculate_diffs() function. If you want to permute by a different variable you need to modify the function.

### Do the calculation on the  observed values
observed_diffs.wang.wangtype.rA.allsamples<-raw_relatedness_allsamples_wang_bothways.sametide %>% 
  dplyr::filter(wangtype == "rA") %>% 
  dplyr::select(c(Tide,Wang)) %>% 
  calculate_diffs()

## map our function over all elements in the permutation list that was output by permute() to test the differences between the means in Tidal versus NonTidal environments on all the permuted values - this step can be time consuming if you have lots of permutations

permuted_diffs_list.wangtype.rA.allsamples <- map(permutations.wangtype.rA.allsamples$perm, calculate_diffs) 

## Pull differences from permuted data out of the list and into a data_frame
permuted_diffs_df.wangtype.rA.allsamples <-do.call("rbind", permuted_diffs_list.wangtype.rA.allsamples) 

# map our calculate_p_value() function over the distribution of permuted differences to calculate the number of permutations that have a more extreme difference than our observed value. This number out of the number of permutations is the p-value.

p_values_tibble.wang.wangtype.rA.allsamples <- 
  map2(names(observed_diffs.wang.wangtype.rA.allsamples), 
       observed_diffs.wang.wangtype.rA.allsamples, 
       ~ {calculate_p_value(.x, .y, 
                            permuted_diffs_df.wangtype.rA.allsamples[[.x]])
                                 }) %>%
  bind_rows()

#++++++++++++++++++++++++++++++++++++++++
#### Permute for rW for all samples.   #### 
#++++++++++++++++++++++++++++++++++++++++

## set up the permutation list.

permutations.wangtype.rW.allsamples <- raw_relatedness_allsamples_wang_bothways.sametide %>% 
  dplyr::filter(wangtype == "rW") %>% 
  modelr::permute(., 10000, Wang)

## Use our calculate_diffs() function to calculate differences between means of values in NonTidal versus Tidal environments.  The variable Tide is hard-wired into the calculate_diffs() function. If you want to permute by a different variable you need to modify the function.

### Do the calculation on the  observed values
observed_diffs.wang.wangtype.rW.allsamples<-
  raw_relatedness_allsamples_wang_bothways.sametide%>% 
  dplyr::filter(wangtype == "rW") %>% 
  dplyr::select(c(Tide,Wang)) %>% 
  calculate_diffs()

## map our function over all elements in the list that was output by permute() to test the differences between the means in Tidal versus NonTidal environments on all the permuted values - this step can be time consuming if you have lots of permutations

permuted_diffs_list.wangtype.rW.allsamples <- map(permutations.wangtype.rW.allsamples$perm, calculate_diffs) 

## Pull differences from permuted data out of the list and into a data_frWme
permuted_diffs_df.wangtype.rW.allsamples <-do.call("rbind", 
                                                   permuted_diffs_list.wangtype.rW.allsamples) 

# map our calculate_p_value() function over the distribution of permuted differences to calculate the number of permutations that have a more extreme difference than our observed value. This number out of the number of permutations is the p-value.

p_values_tibble.wang.wangtype.rW.allsamples <- 
  map2(names(observed_diffs.wang.wangtype.rW.allsamples), 
       observed_diffs.wang.wangtype.rW.allsamples, 
       ~ {calculate_p_value(.x, .y, permuted_diffs_df.wangtype.rW.allsamples[[.x]])
                                                    }) %>%
  bind_rows()


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Correlations with River Km and Within Site Means of rW and rA ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++
#### Read in RiverKm ####
#+++++++++++++++++++++

## Read data file of river kilometers needed for correlations.

RiverKm <- readr::read_csv(here("data","PR.Pop.Centroids.UTMS.csv")) %>% 
  dplyr:: select(NewPop, RiverKm)

#++++++++++++++++++++++++++++++++++++++++++++++++++++
#### Read in Wang data and join river kilometers ----
#++++++++++++++++++++++++++++++++++++++++++++++++++++

wang.allpotomac.within.and.among.pops <- readr::read_csv(here("processed.data",
                                                              "related", 
                                                              "allsamples",
                                                              "wang.rW.and.rA.allsamples.by.site.csv")) %>% 
  dplyr::left_join(.,RiverKm, by= "NewPop")

## Choose response variables

vars.wang <-c("rW.mean", "rA.mean")

#+++++++++++++++++++++++++++++++
#### Full River Correlation ####
#+++++++++++++++++++++++++++++++

corr.results.fullriver.wang.allsample <- map(vars.wang, 
                                  corrfunct, 
                                  df=wang.allpotomac.within.and.among.pops) %>% 
  purrr::set_names(.,vars.wang) %>% 
  purrr::map_dfr(broom::tidy) %>% 
  dplyr::select(-column1, -n) %>% 
  dplyr::rename(r.full.river=estimate, p.full.river =p.value)

#++++++++++++++++++++++++++++++++
#### Tidal Only Correlation ####
#+++++++++++++++++++++++++++++++

corr.results.Tidal.wang.allsample <- map(vars.wang, 
                              corrfunct, 
                              df=filter(wang.allpotomac.within.and.among.pops, 
                                        Tide=="Tidal")) %>% 
  purrr::set_names(.,vars.wang) %>% 
  purrr::map_dfr(broom::tidy) %>% 
  dplyr::select(-column1, -n) %>% 
  dplyr::rename(r.tidal=estimate, p.tidal =p.value)

#+++++++++++++++++++++++++++++++++++
####  NonTidal Only Correlation ####
#+++++++++++++++++++++++++++++++++++

corr.results.Non.Tidal.wang.allsample <- map(vars.wang, 
                                  corrfunct, 
                                  df=filter(wang.allpotomac.within.and.among.pops, 
                                            Tide=="NonTidal")) %>% 
  purrr::set_names(.,vars.wang) %>% 
  purrr::map_dfr(broom::tidy) %>% 
  dplyr::select(-column1, -n) %>% 
  dplyr::rename(r.nontidal=estimate, p.nontidal =p.value)


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#        UNIQUE MLGs      ####
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#+++++++++++++++++++
#### Read in data ####
#+++++++++++++++++++

#  Read in and process the raw relatedness pairwise MLG results.Includes rw and rA for sample-level comparisons. Only MLGs from the same tidal regime are retained.

raw_relatedness_noreps_wang_bothways.sametide <- 
  readr::read_csv(here("processed.data",
                       "related","noreps", 
                       "raw_relatedness_noreps_wang_bothways.csv")) %>% 
  filter(Tide.Comp == "Same.Tide") 


#++++++++++++++++++++
## Permutations ####
#++++++++++++++++++++

#++++++++++++++++++++++++++++++++++++++++++++++++++++++
##### Permute Overall r between Tidal and NonTidal  ----
#++++++++++++++++++++++++++++++++++++++++++++++++++++++

## permute creates the chosen number of dataframes in each of which each chosen variable has been permuted. This action creates a list in which each element is a permuted data frame.

permutations.wang.noreps<- modelr::permute(raw_relatedness_noreps_wang_bothways.sametide, 10000, Wang)

## Use our calculate_diffs() function to calculate differences between means of values in NonTidal versus Tidal environments.  The variable Tide is hard-wired into the calculate_diffs() function. If you want to permute by a different variable you need to modify the function.

### Do the calculation on the  observed values
observed_diffs.wang.noreps<-raw_relatedness_noreps_wang_bothways.sametide %>% 
  dplyr::select(c(Tide,Wang)) %>% 
  calculate_diffs()

## map our function over all elements in the list that was output by permute() to test the differences between the means in Tidal versus NonTidal environments on all the permuted values - this step can be time consuming if you have lots of permutations

permuted_diffs_list.wang.noreps <- map(permutations.wang.noreps$perm, calculate_diffs) 

## Pull differences from permuted data out of the list and into a data_frame
permuted_diffs_df.wang.noreps <-do.call("rbind", permuted_diffs_list.wang.noreps) 

# map our calculate_p_value() function over the distribution of permuted differences to calculate the number of permutations that have a more extreme difference than our observed value. this number out of the number of permutations is the p-value.

p_values_tibble.wang.noreps <- 
  map2(names(observed_diffs.wang.noreps), 
       observed_diffs.wang.noreps, 
       ~ {calculate_p_value(.x, .y, 
                            permuted_diffs_df.wang.noreps[[.x]])
       }) %>%
  bind_rows()



#++++++++++++++++++++++++++++++++++++++++++++++++
##### Permute rA for all MLGs.   #### 
#++++++++++++++++++++++++++++++++++++++++++++++++

permutations.wangtype.rA.noreps <- raw_relatedness_noreps_wang_bothways.sametide %>% 
  dplyr::filter(wangtype == "rA") %>% 
  modelr::permute(., 10000, Wang)

## Use our calculate_diffs() function to calculate differences between means of values in NonTidal versus Tidal environments.  The variable Tide is hard-wired into the calculate_diffs() function. If you want to permute by a different variable you need to modify the function.

### Do the calculation on the  observed values
observed_diffs.wang.wangtype.rA.noreps<-
  raw_relatedness_noreps_wang_bothways.sametide %>% 
  dplyr::filter(wangtype == "rA") %>% 
  dplyr::select(c(Tide,Wang)) %>% 
  calculate_diffs()

## map our function over all elements in the list that was output by permute() to test the differences between the means in Tidal versus NonTidal environments on all the permuted values - this step can be time consuming if you have lots of permutations

permuted_diffs_list.wangtype.rA.noreps <- 
  map(permutations.wangtype.rA.noreps$perm, calculate_diffs) 

## Pull differences from permuted data out of the list and into a data_frame
permuted_diffs_df.wangtype.rA.noreps <-do.call("rbind", permuted_diffs_list.wangtype.rA.noreps) 

# map our calculate_p_value() function over the distribution of permuted differences to calculate the number of permutations that have a more extreme difference than our observed value. this number out of the number of permutations is the p-value.

p_values_tibble.wang.wangtype.rA.noreps <- 
  map2(names(observed_diffs.wang.wangtype.rA.noreps), 
       observed_diffs.wang.wangtype.rA.noreps, 
       ~ {calculate_p_value(.x, .y, 
                            permuted_diffs_df.wangtype.rA.noreps[[.x]])
       }) %>%
  bind_rows()

#++++++++++++++++++++++++++++++++++++
##### Permute rW for all MLGs.   #### 
#+++++++++++++++++++++++++++++++++++++

permutations.wangtype.rW.noreps <- raw_relatedness_noreps_wang_bothways.sametide %>% 
  dplyr::filter(wangtype == "rW") %>% 
  modelr::permute(., 10000, Wang)

## Use our calculate_diffs() function to calculate differences between means of values in NonTidal versus Tidal environments.  The variable Tide is hard-wired into the calculate_diffs() function. If you want to permute by a different variable you need to modify the function.

### Do the calculation on the  observed values
observed_diffs.wang.wangtype.rW.noreps<-raw_relatedness_noreps_wang_bothways.sametide %>% 
  dplyr::filter(wangtype == "rW") %>% 
  dplyr::select(c(Tide,Wang)) %>% 
  calculate_diffs()

## map our function over all elements in the list that was output by permute() to test the differences between the means in Tidal versus NonTidal environments on all the permuted values - this step can be time consuming if you have lots of permutations

permuted_diffs_list.wangtype.rW.noreps <- map(permutations.wangtype.rW.noreps$perm, calculate_diffs) 

## Pull differences from permuted data out of the list and into a data_frWme
permuted_diffs_df.wangtype.rW.noreps <-do.call("rbind", permuted_diffs_list.wangtype.rW.noreps) 

# map our calculate_p_value() function over the distribution of permuted differences to calculate the number of permutations that have a more extreme difference than our observed value. this number out of the number of permutations is the p-value.

p_values_tibble.wang.wangtype.rW.noreps <- 
  map2(names(observed_diffs.wang.wangtype.rW.noreps), 
       observed_diffs.wang.wangtype.rW.noreps, 
       ~ {calculate_p_value(.x, .y, permuted_diffs_df.wangtype.rW.noreps[[.x]])
       }) %>%
  bind_rows()

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Correlations with River Km and Within Sites Means of rW and rA ####
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++
## Read in data ####
#+++++++++++++++++++

#The river km data should have been read in already - check to see that it has and read it in if the object does not exist.

if (!exists("RiverKm")) {
  # If it doesn't exist, read the data file
  RiverKm <- read_csv(here("data", "RiverKm.csv"))
}

wang.noreps.within.and.among.pops<- 
  readr::read_csv(here("processed.data",
                       "related",
                       "noreps",
                       "wang.rW.and.rA.noreps.relatedness.by.site.csv")) %>% 
  dplyr::left_join(.,RiverKm, by= "NewPop")


## Choose response variables
vars.wang <-c("rW.mean", "rA.mean")

#+++++++++++++++++++++++++++++++++
####  Full River Correlation ####
#+++++++++++++++++++++++++++++++++

corr.results.fullriver.wang.noreps <- purrr::map(vars.wang, 
                                             corrfunct, 
                                             df=wang.noreps.within.and.among.pops) %>% 
  purrr::set_names(.,vars.wang) %>% 
  purrr::map_dfr(broom::tidy) %>% 
  dplyr::select(-column1, -n) %>% 
  dplyr::rename(r.full.river=estimate, p.full.river =p.value)

#+++++++++++++++++++++++++++++++++++
####  Tidal Only Correlation ####
#+++++++++++++++++++++++++++++++++++

corr.results.Tidal.wang.noreps <- 
  map(vars.wang, 
      corrfunct, 
      df=filter(wang.noreps.within.and.among.pops, 
                Tide=="Tidal")) %>% 
  purrr::set_names(.,vars.wang) %>% 
  purrr::map_dfr(broom::tidy) %>% 
  dplyr::select(-column1, -n) %>% 
  dplyr::rename(r.tidal=estimate, p.tidal =p.value)

#+++++++++++++++++++++++++++++++++++
####  NonTidal Only Correlation ####
#+++++++++++++++++++++++++++++++++++

corr.results.Non.Tidal.wang.noreps <- 
  map(vars.wang, 
      corrfunct, 
      df=filter(wang.noreps.within.and.among.pops, 
                Tide=="NonTidal")) %>% 
  purrr::set_names(.,vars.wang) %>% 
  purrr:: map_dfr(broom::tidy) %>% 
  dplyr::select(-column1, -n) %>% 
  dplyr::rename(r.nontidal=estimate, p.nontidal =p.value)
