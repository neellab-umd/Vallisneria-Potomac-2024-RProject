## Author: Maile Neel
## 2024

## This code generates genind and genpop objects for the Potomac using the data on microsatellite genotypes and site information in the file .data/allpotomac.microsatellite.data.csv. These objects are used in analyses of population genetic diversity. The code is designed to work within the Potomac.Manuscript.2024.R project.

## The different population analyses are done in separate R scripts in this project.

## The objects created by this script need to exist for the other scripts to work. If these three rds files exist in the ./processed.data folder you do not need to run this script again. 
## list.of.genind.objects.rds
## list.of.genpop.objects.rds.

## The related package takes different data formats; data for those analyses are read in via code the relevant scripts.

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
## Install & load packages and Source Functions  ####
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,
               here,
               purrr,
               adegenet,
               update = FALSE)


## Load our custom function from the functions folder in this project - first check to see if it is loaded, and if not, load it.

if(!exists("create_geninds", mode="function")) source(here::here("functions", "create_geninds.r"))

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
## Load Microsatellite Data ######
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

# Read data file in as a dataframe - will be
## converted to a genind object later.

gendata.allpotomac<-read_csv(here("data","allpotomac.microsatellite.data.csv")) %>% 
  dplyr::select(-c(Sample_related,TideCode,Collector,aagx030.1:m16.2)) %>% 
  drop_na(Clone.ID.2018)

## check that the data file was read in correctly.
dim(gendata.allpotomac)
names(gendata.allpotomac)

## check subset of the content
## display the first 9 columns for 50 rows and the last 6 rows
gendata.allpotomac[200:250,1:9]
tail(gendata.allpotomac[,16:20])

### Only NonTidal Sites ----
gendata.allpotomac.NonTidal<-gendata.allpotomac %>% 
  dplyr::filter(Tide=="NonTidal")
dim(gendata.allpotomac.NonTidal)

### Only Tidal Sites ----
gendata.allpotomac.Tidal<-gendata.allpotomac %>% 
  filter(Tide=="Tidal")
dim(gendata.allpotomac.Tidal)

## create a dataframe with only one instance of each unique MLG within each site in the base dataframe.

gendata.noreps.potomac<-gendata.allpotomac %>% 
  dplyr::distinct(NewPop, Clone.ID.2018, .keep_all=TRUE) 

dim(gendata.noreps.potomac)

## Only Tidal Sites
gendata.noreps.potomac.Tidal<-gendata.noreps.potomac %>% 
  dplyr::filter(Tide=="Tidal")
dim(gendata.noreps.potomac.Tidal)

## Only NonTidal Sites
gendata.noreps.potomac.NonTidal<-gendata.noreps.potomac %>% 
  dplyr::filter(Tide=="NonTidal")
dim(gendata.noreps.potomac.Tidal)

## Only one instance of each MLG overall - for QA/QC testing

gendata.noreps.atall<-gendata.allpotomac %>% 
  dplyr::distinct(Clone.ID.2018, .keep_all=TRUE) 

### Create one list to hold all the gendata dataframes ----
list.of.gendata.dataframes<-setNames(map(ls(pattern="^gendata."), get),
                                     ls(pattern="^gendata.")) 
names(list.of.gendata.dataframes)


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
##  Create adegenet genind Objects ####
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

## Convert all our dataframes to genind objects using our custom function (create_geninds) that uses df2genind. You are using a loop to create genind.objects for each of the raw data sets you need to summarize.

## Map our custom function over the list of gendata objects (make sure you sourced the function when you loaded packages).

list.of.genind.objects <- purrr::map(list.of.gendata.dataframes, ~create_geninds(.x))

## change the names to include genind instead of gendata.

names(list.of.genind.objects)<- sub("gendata.", "genind.", names(list.of.genind.objects))

## Make a couple of other genind objects using Tide as the population and the strata.

## creates an all_reps ind object that uses Tidal/NonTidal as the population labels

genind.allpotomac.bytide<- list.of.genind.objects$genind.allpotomac

pop(genind.allpotomac.bytide) <- genind.allpotomac.bytide@other$Tide

adegenet::strata(genind.allpotomac.bytide) <- 
  data.frame(adegenet::pop(genind.allpotomac.bytide))

adegenet::nameStrata(genind.allpotomac.bytide) <- ~Tide

list.of.genind.objects$genind.allpotomac.bytide <- genind.allpotomac.bytide

## createa noreps ind object that uses Tidal/NonTidal as the population labels and sets it as the strata as well.

genind.noreps.potomac.bytide<-list.of.genind.objects$genind.noreps.potomac

adegenet::pop(genind.noreps.potomac.bytide)<-genind.noreps.potomac.bytide@other$Tide

adegenet::strata(genind.noreps.potomac.bytide)<-
  data.frame(adegenet::pop(genind.noreps.potomac.bytide))
adegenet::nameStrata(genind.noreps.potomac.bytide) <- ~Tide

list.of.genind.objects$genind.noreps.potomac.bytide <- genind.noreps.potomac.bytide

saveRDS(list.of.genind.objects, file = here::here("processed.data","list.of.genind.objects.rds"))

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
## Create adegenet genpop Objects ####
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

## Convert each of the genind objects into adegenet genpop objects. process.other = TRUE keeps all the @other slots with the populations.

## To get map to process the @other slots, I created the following function for genind2genpop that buries the internal summarizing into a single argument function.

map.convert.to.genpop <- function(x) {
  adegenet::genind2genpop(x, process.other = TRUE, other.action = mean)
}

## map that function over all the genind objects in your list
list.of.genpop.objects <-purrr::map(list.of.genind.objects, ~map.convert.to.genpop(.x))

## change the names to include genpop instead of genind.
names(list.of.genpop.objects)<- sub("genind.", "genpop.", names(list.of.genpop.objects))

saveRDS(list.of.genpop.objects, file = here::here("processed.data","list.of.genpop.objects.rds"))

## All your dataframes, genind, and genpop objects and have now been created. They are processed  using other scripts in the project using the packages and adegenet, hierfstat, and MMOD.

## Delete the dataframes with names that start with "gendata" from the environment to tidy up.

#rm(list=ls(pattern="^gendata."))
