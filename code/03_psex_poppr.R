# Maile Neel 2024

## Usage ----

## This script contains the poppr code used run the pgen and pex analysis.  It is designed to be run within the R Project 01_Potomac.Manuscript.RProject.Rproj.

#It uses the genind objects created in the 01_PR_Manuscript_adegenet_Create_Dataframes.R script. That script will save the list of geninds to an rds file called list.of.genind.objects.rds in the processed.data folder. If that file exists, you do not need to run that script again and you are ready to run this one.

## If not, you will need to run 01_PR_Manuscript_adegenet_Create_Dataframes.R. 

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Install and Load Packages ####
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

if (!require("pacman")) install.packages("pacman")
# installs missing packages that are found on CRAN, and then loads the packages
pacman::p_load(here,
               tidyverse,
               magrittr,
               poppr
               )
#XXXXXXXXXXXXXXXXXXXXXXXXXXX
# Load custom functions ####
#XXXXXXXXXXXXXXXXXXXXXXXXXXX

## This script requires custom functions that are in the functions subfolder of the project.

if(!exists("create.genclones", mode="function")) source(here::here("functions", "create.genclones.r"))

if(!exists("pgen.and.sex", mode="function")) source(here::here("functions", "create.psex.and.pgen.summary.r"))

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Read in lists of genind objects ####
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

## These lists were created created in the script that includes in its name the string *_adegenet_create_dataframes.r

list.of.genind.objects <- 
  list.files("./processed.data/", 
             pattern=".genind.objects.rds", 
             all.files=TRUE,
             full.names=TRUE) %>% 
  readRDS(.) 

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
## Process genind objects using poppr ####
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

## Code to convert all the genind objects in your list of genclone objects that are needed by poppr specify using my mlls coming from the CLONE.ID.2018 slot in the genind object to keep it from calling new mlls

## Convert only the full dataset geninds to a genclone object by negating the string "noreps" and the within CRV data set.

list.of.genclone.objects <- names(list.of.genind.objects) %>% 
  stringr::str_detect('noreps', negate = TRUE) %>%
  purrr::keep(list.of.genind.objects, .) %>% 
  purrr::map(., ~create.genclones(.))

names(list.of.genclone.objects)<- sub("genind.", "genclone.", names(list.of.genclone.objects))

names(list.of.genclone.objects)

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Checking Ability to Detect Sexual Reproduction ####
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

# Here we explore the ability to detect if mlgs arose from sexual reproduction versus clonal reproduction given the allelic diversity in our microsatellites.

### Plot power of loci ####
# First plot the number of mlgs detected with increasing numbers of loci, for all regions combined, and then separately for PR, NY, and ME.

list.of.genclone.objects$genclone.allpotomac %>% 
  poppr::genotype_curve(
    sample = 100,
  maxloci = 0,
  quiet = FALSE,
  thresh = 1,
  plot = TRUE,
  drop = TRUE,
  dropna = TRUE
)

#+++++++++++++++++++++++++++++++++++++++++
######## Calculate pgen and psex #########
#+++++++++++++++++++++++++++++++++++++++++

## Use the custom function loaded above to run poppr functions pgen(), psex() for method = "single" and method ="multiple" on a data set, gather up the results into a dataframe, and output the results to a csv file.  Again, make sure you have sourced the function if you want to use it.

## pgen() is the probability of encountering a genotype by chance given the allele frequencies.
## psex() gives the probability of encountering a genotype more than once by chance.
## If psex() is calculated using method = "single", the result is the probability of seeing the mlg a second time due to sexual reproduction. 

## If method = "multiple" the result is the probability of sexual reproduction yielding the mlg the number of times it was seen. 

## Details of the calculations are in the poppr manual.

## We do the calculations for the overall river and by separate populations. Both were used to evaluate probabilities of multiple occurrences of MLGs.  If you calculate psex by population, the first instance of an MLG in that population will have a probability near 1 even though pgen is way smaller than that.  If you calculate without considering population, the probability of multiple occurrences is calculated correctly, but I am not sure from the poppr manual that Fis is taken into account.

## The function uses the default correction for rare alleles (i.e., correcting to 1/n). That could be changed by modifying the function in the file create.psex.and.pgen.summary.r.

## Run the function on one genclone object ####

allpotomac.pgen <- pgen.and.sex(list.of.genclone.objects$genclone.allpotomac)

# You can also run the function on a list of genclone objects ####

## Summarize Results ----

# psex method single for mlgs with two observations
Num.Psex.gt.0.05 <- allpotomac.pgen %>% 
  dplyr::filter(n == 2) %>% 
  dplyr::group_by(pop) %>% 
  dplyr::summarize(Single = length(mlg[Psex.second >= .01]))

# psex method multiple for the number of times each mlg was seen - include all mlgs seen 2 or more times.
Num.Pmult.gt.0.05 <- allpotomac.pgen  %>% 
  filter(n >= 2) %>% 
  #the first instance of each mlg has a prob of ~0.99.Filter those out and make sure the rest are < 0.05
  dplyr::group_by(mlg) %>% 
  dplyr::filter(Psex.multiple < 0.95) %>% 
  dplyr::summarize(multiple = length(mlg[Psex.multiple >= .01]))

