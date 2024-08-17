## Maile Neel
## 2024

## Usage ----

## This script contains the poppr code used to generate genotypic diversity statistics. It is designed to be run within the RProject called Potomac.Manuscript.2024.RProject.Rproj.

## It uses the list of genind objects created in the 01.2018_PR_Manuscript_adegenet_create_dataframes.R script. If that script has been run, you will find an rds file called list.of.genind.objects.rds in the processed.data folder.  You can read that in (code is below) to get the genind objects you need.

## If the list.of.genind.objects.rds does not exist, 01.2018_PR_Manuscript_adegenet_create_dataframes.R must be run first.

## The end result of this script is two dataframes. One dataframe has all the MLG information for each tidal environment and the whole river, and one dataframe has the same information from the site level.  These dataframes are saved to the ./processed.data folder as poppr.tide.level.output.allpotomac.csv and poppr.site.level.output.allpotomac.csv. If those csv files are present, this script does not need to be rerun and you can proceed right to the analysis script (02b_Summarize_and_Analyze_poppr_Output.R).

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
## Load Packages & Source Functions ####
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,
               here,
               poppr,
               update = FALSE)

## Source custom functions needed to create genclone objects and to summarize poppr output.

if(!exists("create.genclones", mode = "function"))
  source(here::here("functions", "create.genclones.r"))

if (!exists("create.poppr.tables", mode = "function"))
  source(here::here("functions", "create.poppr.tables.r"))

if (!exists("overall.poppr.table", mode = "function"))
  source(here::here("functions", "overall.poppr.table.r"))


## poppr needs data to be in genclone format. This script will create genclone objects from the genind objects for all samples that were created in the script 01_Create_adegenet_geninds.R using our custom create.genclones function. If there is an rds file called list.of.genclone.objects.rds in the ./processed.data folder, you do not need to create it again. If that rds file is not present, you must use this script to create the genclone objects before proceeding to the Calculate Statistics section of this scrip.

list.of.genclone.objects<-readRDS(here("processed.data", "list.of.genclone.objects.rds")
)

## The function create.genclones() hardwires using MLLs coming from the CLONE.ID.2018 slot in the genind object.  To use with a different dataset in which mlg ids are in a different column, the function needs to be modified.

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Read in genind objects ----
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXX

## Read in rds file with genind objects. the list.of.genind.objects.rds file is not in the processed.data. folder, you need to create it with the script 01_Create_adegenet_geninds.R.

list.of.genind.objects <-
  readRDS(here("processed.data", "list.of.genind.objects.rds"))

## From the list.of.genind.objects you use only the genind object for all samples (the noreps sets are not relevant for MLG diversity analyses).

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Create genclone objects ----
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXX

## map our custom create.genclones() function over only the genind with all samples - denoted in their names as 'allpotomac'. 

list.of.genclone.objects <- names(list.of.genind.objects) %>%
  str_detect('allpotomac') %>%
  purrr::keep(list.of.genind.objects, .) %>%
  map(., ~ create.genclones(.x))

## change the names of items in the output list to include genclone instead of genind.

names(list.of.genclone.objects) <-
  sub("genind.", "genclone.", names(list.of.genclone.objects))

## Write out the list of genclone objects for future use.

write_rds(list.of.genclone.objects,
          here("processed.data", "list.of.genclone.objects.rds"))


#XXXXXXXXXXXXXXXXXXXXXXXXXX
# Calculate Statistics ----
#XXXXXXXXXXXXXXXXXXXXXXXXXX

## Use our custom function to run poppr on the genclone.allpotomac gc object and create summary numbers and the values for table 1. Even if you are not rarefying, you need to give a number for sample size in the data set (n.rare=).

### By Tidal Regime ----

## First run it using Tide as the group to get information on MLGs independent of sites.

overall.mlg.diversity <-
  overall.poppr.table(
    list.of.genclone.objects$genclone.allpotomac,
    group = "Tide",
    rarefy = FALSE,
    n.rare = 5
  )

### By Site ----

## Use NewPop as the group to get population level mlg information. 

## The poppr summary table only rarefies the estimate of eMLG.  To get rarefied values for the other statistics you need to use this code to sample the specified number of samples from each population. Even if you are not rarefying, you need to give the smallest sample size in the data set so our custom function requires it as an argument. 

mlg.diversity.by.site.rarefied <-
  create.poppr.tables(list.of.genclone.objects$genclone.allpotomac,
                      group = "NewPop",
                      rarefy = TRUE,
                      n.rare = 5
  )

## The create.poppr.tables function cretes the variables we use from the original poppr output.  In that output G is effective Simpsons, H is Shannon's. We convert H to the effective version by taking exp(H) and we convert G to its effective form by taking 1/G   eMLG is the number of expected MLG at the smallest sample size 5 based on rarefaction.

## To keep the OrderPop variable with the corresponding NewPop, we need to pull those values out of the @other slot from the original genind and get the unique values.

mlg.diversity.by.site.with.labels <-  as.data.frame(list.of.genind.objects$genind.allpotomac@strata) %>% 
  dplyr::group_by(NewPop) %>% 
  dplyr::summarize(OrderPop = first(OrderPop),
                   Tide = first(Tide)) %>% 
  dplyr::left_join(mlg.diversity.by.site.rarefied, by = "NewPop")

#XXXXXXXXXXXXXXXXXXXXXXXX
# Save poppr Output ----
#XXXXXXXXXXXXXXXXXXXXXXXX

## Write the poppr output dataframe to csv file that can be used to pull summary numbers by tidal regime. 

readr::write_csv(overall.mlg.diversity,
                 here::here(
                   "processed.data",
                   "poppr.tide.level.output.allpotomac.csv"))

## Write the poppr output dataframe to csv file that can be used for summarizing and testing for differences between based on sites in the 02b_Summarize_and_Analyze_poppr_Output.R script. 

readr::write_csv(mlg.diversity.by.site.with.labels,
                                  here::here(
                                    "processed.data",
                                    "poppr.site.level.output.allpotomac.csv"))

### Analysis of the poppr output is done in the script titled  02b_Summarize_and_Analyze_poppr_Output.R  With the mlg.diversity.by.site.with.labels in your environment you can proceed right to that script.  That script has code to read the csv file back in so you do not have to rerun this whole script. 
