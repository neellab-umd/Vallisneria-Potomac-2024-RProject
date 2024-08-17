## Author: Maile Neel
## 2024

## This script is designed to work in the Potomac.Manuscript.2024.R.Project.RProj.  

## It creates the summary dataframes of Wang's r created using the related package as implemented in script 06a_Calculate_Relatedness.R.  If coancestry.results.allsamples.RDATA exists in the ./processed.data/related/allsamples folder and coancestry.results.noreps.RDATA is in the ./processed.data/related/ no reps folder you do not need run that script again. 

## The overall and mean values of Wang's R and values of Wang'r within versus among sites presented in the manuscript and in Figure 4 are generated in this script.  

## Significance testing through permutations of tidal versus non-tidal environments and correlation with river kilometer are implemented in 06c_Relatedness_Permutation.R

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# INSTALL AND LOAD PACKAGES ####
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyr,
               tidyverse,
               here)

#XXXXXXXXXXXXXXXX
# READ DATA ####
#XXXXXXXXXXXXXXXX

## Read in the microsatellite datafile to get the OrderPop and tide information for each ste that will be used for rejoining OrderPop and Tide back after the relatedness data are processed.

potomac_for_joining <- read_csv(here("data","allpotomac.microsatellite.data.csv")) %>% 
  #  drop_na(Clone.ID.2018) %>% 
  dplyr::select(c(IDName,NewPop,OrderPop, Tide, TideCode)) %>% 
  dplyr::group_by(NewPop) %>% 
  dplyr::summarize(OrderPop = first(OrderPop),
                   TideCode = first(TideCode),
                   Tide = first(Tide))

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#            # ALL SAMPLE COANCESTRY SUMMARY       ###########
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

## Load the coancestry.results.all.samples list that was output by coancestry() in the script 06a.Calculate_Relatedness.R

load(here("processed.data","related","allsamples","coancestry.results.allsamples.RDATA"))

##Get the individual labels and the relatedness estimator from the coancestry output object, add column names, and make Wang's r numeric.

raw_relatedness_allsamples_wang <-
  as.data.frame(cbind(coancestry.results.allsamples$relatedness$ind1.id,
                      coancestry.results.allsamples$relatedness$ind2.id,
                      coancestry.results.allsamples$relatedness$wang)) %>% 
  purrr::set_names(c("From","To","Wang")) %>% 
  dplyr::mutate(Wang = as.numeric(Wang))

#Pick up both directions of the comparisons - all froms and all to's. Reverse the order of the To and From columns, but relabel them From and To

raw_relatedness_allsamples_wang_reverse <-raw_relatedness_allsamples_wang %>% 
  relocate(To, From, Wang) %>% 
  purrr::set_names(c("From","To","Wang")) 

## Bind both directions together and split out the components of the From and To labels. Also create Within and Among Site and Tide variables to make summarizing by group easier.

raw_relatedness_allsamples_wang_bothways<-rbind(raw_relatedness_allsamples_wang,
                                        raw_relatedness_allsamples_wang_reverse) %>% 
  tidyr::separate_wider_delim(To, ".", names = c("ToTideCode","ToSite","ToSamp")) %>% 
  tidyr::separate_wider_delim(From, ".", names = c("FromTideCode","FromSite","FromSamp")) %>% 
  dplyr::mutate(Site.Comp = case_when(ToSite == FromSite ~ "Same.Site",
                                        TRUE ~ "Different.Site"),
                Tide.Comp = case_when(ToTideCode == FromTideCode ~ "Same.Tide",
                                      TRUE ~ "Different.Tide")) %>% 
  dplyr::mutate(Tide = case_when(FromTideCode == "AA" ~"NonTidal",
                                 FromTideCode == "AB" ~"Tidal")) %>% 
  dplyr::mutate(wangtype = case_when(Site.Comp == "Same.Site" ~"rW",
                                     Site.Comp == "Different.Site" ~"rA"))

dim(raw_relatedness_allsamples_wang_bothways)

#you now have every comparison in twice, but you only summarize from one direction so you will summarize all comparisons for each individual correctly.

#++++++++++++++++++++++++++++++++++++++++++++++++++++
## Save pairwise r for later use in permutations ----
#+++++++++++++++++++++++++++++++++++++++++++++++++++

write_csv(raw_relatedness_allsamples_wang_bothways, 
          here("processed.data","related","allsamples",
               "raw_relatedness_allsamples_wang_bothways.csv"))

#++++++++++++++++++++++++++++++++++++++++++++++
###  Summarize Wangs r Overall and by Tide ####
#++++++++++++++++++++++++++++++++++++++++++++++

## Summarize overall mean Wangs r and sd for the whole river. Use the original dataframe because you are not choosing by only one From or To category so you are not eliminating duplicates.

raw_relatedness_allsamples_wang %>% 
  dplyr::summarise(mean=mean(Wang),sd(Wang),n())

## Summarize overall mean Wang's r and sd for samples from each tidal regime to all other samples in the river. 

raw_relatedness_allsamples_wang_bothways %>% 
  dplyr::group_by(FromTideCode) %>% 
  dplyr::summarise(mean=mean(Wang),sd(Wang),n())

## Summarize overall mean Wang's r and sd for samples from each tidal regime to all other samples the river separated by whether those individuals were in the same or different tidal regime.

raw_relatedness_allsamples_wang_bothways %>% 
  dplyr::group_by(FromTideCode, Tide.Comp) %>% 
  dplyr::summarise(mean=mean(Wang),sd(Wang),n())


#++++++++++++++++++++++++++++++++++++++++
###   Summarize rW and rA by Sample  ####
#++++++++++++++++++++++++++++++++++++++++

## For each sample, calculate the mean pairwise Wang’s r to all other samples within the same site (rW) and to samples from other sites within the same tidal regime (rA)

rW.and.rA.allsamples.by.sample <- raw_relatedness_allsamples_wang_bothways %>% 
  dplyr:: filter(Tide.Comp != "Different.Tide") %>% ## keep only samples from the same tide
  dplyr::group_by(FromSamp, Site.Comp) %>% 
  dplyr::summarise(mean.Wang = mean(Wang),
                   FromTideCode = first(FromTideCode),
                   FromSite = first(FromSite)) %>% 
  tidyr::pivot_wider(id_cols = c(FromSamp,FromSite,FromTideCode), 
              names_from = Site.Comp, 
              values_from = c(mean.Wang)) %>% 
  dplyr::rename(Wang.allreps.rA = Different.Site, Wang.allreps.rW = Same.Site)  %>%
  dplyr::left_join(potomac_for_joining, by = join_by(FromSite == NewPop)) %>% 
  dplyr:: mutate(Tide = case_when(FromTideCode == "AA" ~"NonTidal",
                                  FromTideCode == "AB" ~"Tidal")) %>% 
  dplyr::select(-c(FromTideCode)) 

write_csv(rW.and.rA.allsamples.by.sample,"./processed.data/related/allsamples/wang.rW.and.rA.allsamples.by.sample.csv")

## Calculate the mean rA and rW for samples from the tidal versus nontidal regimes and over the whole river.

rw.and.ra.allsamples.samples.bytide <- rW.and.rA.allsamples.by.sample %>% 
  bind_rows(., rW.and.rA.allsamples.by.sample %>% mutate(FromTideCode = as.character("All"))) %>%
  dplyr:: group_by(FromTideCode) %>% 
  dplyr::summarise(mean.rW = mean(Wang.allreps.rW),sd.rW = sd(Wang.allreps.rW),
                   mean.rA = mean(Wang.allreps.rA),sd.rA = sd(Wang.allreps.rA),
                   n()
                   )

#+++++++++++++++++++++++++++++++++++++++
###   Summarize rW and rA by Site   ####
#+++++++++++++++++++++++++++++++++++++++

rW.and.rA.allsamples.by.site <- rW.and.rA.allsamples.by.sample %>% 
  dplyr::group_by(FromSite) %>% 
  dplyr::summarise(rW.mean=mean(Wang.allreps.rW), 
                   rW.sd=sd(Wang.allreps.rW),
                   rA.mean=mean(Wang.allreps.rA), 
                   rA.sd=sd(Wang.allreps.rA)) %>% 
  dplyr::left_join(potomac_for_joining, by = join_by(FromSite == NewPop)) %>% 
  dplyr::rename(NewPop = FromSite) %>% 
  dplyr:: select(-TideCode)


## Calculate means of site-level rW and rA by tide

rW.and.rA.allsamples.sites.bytide <- rW.and.rA.allsamples.by.site %>% 
  bind_rows(., rW.and.rA.allsamples.by.site %>% mutate(Tide = as.character("All"))) %>%
  dplyr::group_by(Tide) %>% 
  dplyr::summarise(rW = mean(rW.mean),sd.rW = sd(rW.mean),
                   rA = mean(rA.mean),sd.rA = sd(rW.mean),
                   n()
  )

#+++++++++++++++++++++++++++++++++
###   Write rW and rA to csv  ####
#+++++++++++++++++++++++++++++++++

write_csv(rW.and.rA.allsamples.by.site,"./processed.data/related/allsamples/wang.rW.and.rA.allsamples.by.site.csv")

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#          MLG COANCESTRY SUMMARY    ###########
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

## Load the coancestry.results.noreps list that was output by coancestry() in the script 06a.Calculate_Relatedness.R

load(here("processed.data","related","noreps","coancestry.results.noreps.RDATA"))

##Get the individual labels and the relatedness estimator from the coancestry output object, add column names, and make Wang's numeric.

raw_relatedness_noreps_wang <-
  as.data.frame(cbind(coancestry.results.noreps$relatedness$ind1.id,
                      coancestry.results.noreps$relatedness$ind2.id,
                      coancestry.results.noreps$relatedness$wang)) %>% 
  purrr::set_names(c("From","To","Wang")) %>% 
  dplyr::mutate(Wang = as.numeric(Wang))

## Pick up both directions of the comparisons - all froms and all to's. Reverse the order of the To and From columns, but relabel them From and To

raw_relatedness_noreps_wang_reverse <- raw_relatedness_noreps_wang %>% 
  dplyr::relocate(To, From, Wang) %>% 
  purrr::set_names(c("From","To","Wang")) 

## Bind both directions together and split out the components of the From and To labels. Also create Within and Among Site and Tide variables to make summarizing by group easier.

raw_relatedness_noreps_wang_bothways <- rbind(raw_relatedness_noreps_wang,
                                            raw_relatedness_noreps_wang_reverse) %>% 
  tidyr::separate_wider_delim(To, ".", 
                              names = c("ToTideCode","ToSite","ToSamp")) %>% 
  tidyr::separate_wider_delim(From, ".", 
                              names = c("FromTideCode","FromSite","FromSamp")) %>% 
  dplyr::mutate(Site.Comp = case_when(ToSite == FromSite ~ "Same.Site",
                                      TRUE ~ "Different.Site"),
                Tide.Comp = case_when(ToTideCode == FromTideCode ~ "Same.Tide",
                                      TRUE ~ "Different.Tide")) %>% 
  dplyr::mutate(Tide = case_when(FromTideCode == "AA" ~"NonTidal",
                                 FromTideCode == "AB" ~"Tidal")) %>% 
  dplyr::mutate(wangtype = case_when(Site.Comp == "Same.Site" ~"rW",
                                  Site.Comp == "Different.Site" ~"rA"))

dim(raw_relatedness_noreps_wang_bothways)

#you now have every comparison in twice, but you only use it to summarize from one direction so you will summarize all comparisons for each individual correctly.

#+++++++++++++++++++++++++++++++++++++++++++++++
## Save pairwise r for later use in permutations ----
#+++++++++++++++++++++++++++++++++++++++++++++++

write_csv(raw_relatedness_noreps_wang_bothways, 
          here("processed.data","related","noreps",
               "raw_relatedness_noreps_wang_bothways.csv"))

#++++++++++++++++++++++++++++++++++++++++++++++++
###  Summarize Wang's r Overall and by Tide  ####
#++++++++++++++++++++++++++++++++++++++++++++++++

## Summarize overall mean Wangs r and sd for the whole river. Use the original dataframe because you are not choosing by only one From or To category so you are not eliminating duplicates.

raw_relatedness_noreps_wang %>% 
  dplyr::summarise(mean=mean(Wang),sd(Wang),n())

## Summarize overall mean Wang's r and sd for samples from each tidal regime to all other samples in the river. 

raw_relatedness_noreps_wang_bothways %>% 
  dplyr::group_by(FromTideCode) %>% 
  dplyr::summarise(mean=mean(Wang),sd(Wang),n())

## Summarize overall mean Wang's r and sd for samples from each tidal regime to all other samples the river separated by whether those individuals were in the same or different tidal region.

raw_relatedness_noreps_wang_bothways %>% 
  dplyr::group_by(FromTideCode, Tide.Comp) %>% 
  dplyr::summarise(mean=mean(Wang),sd(Wang),n())


#++++++++++++++++++++++++++++++++++++++++
###   Summarize rW and rA by Sample  ####
#++++++++++++++++++++++++++++++++++++++++

## For each sample, calculate the mean pairwise Wang’s r to all other samples within the same site (rW) and to samples from other sites within the same tidal regime (rA)

rW.and.rA.noreps.by.sample <- raw_relatedness_noreps_wang_bothways %>% 
  dplyr::filter(Tide.Comp != "Different.Tide") %>% 
  dplyr::group_by(FromSamp, Site.Comp) %>% 
  dplyr::summarise(mean.Wang = mean(Wang),
                   FromTideCode = first(FromTideCode),
                   FromSite = first(FromSite)) %>% 
  tidyr::pivot_wider(id_cols = c(FromSamp,FromSite,FromTideCode), 
              names_from = Site.Comp, 
              values_from = c(mean.Wang)) %>% 
  dplyr::rename(Wang.noreps.rA = Different.Site, Wang.noreps.rW = Same.Site) %>% 
  dplyr:: mutate(Tide = case_when(FromTideCode == "AA" ~"NonTidal",
                          FromTideCode == "AB" ~"Tidal")) %>% 
  dplyr::select(-c(FromTideCode)) 

write_csv(rW.and.rA.noreps.by.sample,here("processed.data", "related", "noreps", "wang.rW.and.rA.noreps.by.sample.csv"))

## Calculate the mean rA and rW for samples from the tidal versus nontidal regimes and over the whole river.

rw.and.ra.noreps.samples.bytide <- rW.and.rA.noreps.by.sample %>% 
  bind_rows(., rW.and.rA.noreps.by.sample %>% mutate(FromTideCode = as.character("All"))) %>%
  dplyr::group_by(FromTideCode) %>% 
  dplyr::summarise(Wang.noreps.rW.mean = mean(Wang.noreps.rW),sd.rW = sd(Wang.noreps.rW.mean),
                   Wang.noreps.rA.mean = mean(Wang.noreps.rA),sd.rA = sd(Wang.noreps.rA.mean),
                   n()
  )

#+++++++++++++++++++++++++++++++++++++++
###   Summarize rW and rA by Site   ####
#+++++++++++++++++++++++++++++++++++++++

rW.and.rA.noreps.by.site <- rW.and.rA.noreps.by.sample %>% 
  group_by(FromSite) %>% 
  dplyr::summarise(rW.mean = mean(Wang.noreps.rW), 
                   rW.sd = sd(Wang.noreps.rW),
                   rA.mean = mean(Wang.noreps.rA), 
                   rA.sd = sd(Wang.noreps.rA)) %>% 
  dplyr::left_join(potomac_for_joining, by = join_by(FromSite == NewPop)) %>% 
  dplyr:: rename(NewPop = FromSite) %>% 
  dplyr::select(-TideCode)


## Calculate means of site level rW and rA by tide

rW.and.rA.noreps.sites.bytide <- rW.and.rA.noreps.by.site %>% 
  bind_rows(., rW.and.rA.noreps.by.site %>% mutate(Tide = as.character("All"))) %>%
  dplyr::group_by(Tide) %>% 
  dplyr::summarise(rW = mean(rW.mean),sd.rW = sd(rW.mean),
                   rA = mean(rA.mean),sd.rA = sd(rA.mean),
                   n()
  )

#+++++++++++++++++++++++++++++++++
###   Write rW and rA to csv  ####
#+++++++++++++++++++++++++++++++++

write_csv(rW.and.rA.noreps.by.site,"./processed.data/related/noreps/wang.rW.and.rA.noreps.relatedness.by.site.csv")
