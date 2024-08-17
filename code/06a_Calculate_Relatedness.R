## Author: Maile Neel
## 2024

## This script is designed to work in the Potomac.Manuscript.2024.R.Project.RProj.  

## It does all the calculations we use from package related using the allpotomac.microsatellite.data.csv from the data folder as input. 

## The data output from coancestry() are written to files that are processed in the script 06b_Summarize_Relatedness_Output.R. The permutation tests and correlation with river kilometer are done in 06C_Relatedness_Permutation.R

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# INSTALL AND LOAD PACKAGES ----
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,
               devtools,
               here)

#the related package is not on CRAN. You need to get it from github https://github.com/timothyfrasier/related
#Download the tar/zip file, and install it from a tar.gz by selecting "Tools", Install Packages, and Install from Package Archive File.

#You have to have Rtools (matched to your R version) installed to compile the package.
#If you get a message that a make file was not found you don't have Rtools
#You can get that here https://cran.r-project.org/bin/windows/Rtools/

#Note that as of 2024, related is not working properly when R is run under windows. These analyses were run on a linux server.

library(related)

#XXXXXXXXXXXXXXX
# READ DATA ----
#XXXXXXXXXXXXXXX

## Data setup explanation ####

#Related requires a different format for loci and population labels than all the other analyses so you use different columns from the allpotomac.microsatelite.data.csv file and you concatenate a single label for each individual that includes a two-letter group code for Tidal regime, population source, and individual id. related truncates the name so make sure the label doesn't get too long.

## related requires a two letter id code for groups that you want to test.
#We store that code in the variable TideCode in which AA is nontidal and AB is tidal

## I use the IDName column because it does not have periods in the individual IDs. Later you will split the different levels of the label out to columns based on periods so you don't want periods as part of the sample ID.

## All other data will get stripped away but you first will use the Clone.ID.2018 to create a dataset of only unique MLGs within each NewPop.

potomac_for_related <- readr::read_csv(here("data","allpotomac.microsatellite.data.csv")) %>% 
  tidyr::drop_na(Clone.ID.2018) %>% 
  dplyr::select(c(IDName,NewPop,TideCode,aagx030.1:m16.2, Clone.ID.2018)) %>% 
  dplyr::mutate(across(aagx030.1:m16.2, as.numeric),
                across(aagx030.1:m16.2, ~replace_na(.x,0)),
                TidePopSamp = paste0(TideCode,
                               ".", NewPop, 
                               ".", IDName)) %>% 
  dplyr::relocate(TidePopSamp,aagx030.1:m16.2)

#check that data look like they were brought in correctly
#display the first 12 columns for rows 200-250
potomac_for_related[200:250,1:12]
#display a different subset of data
potomac_for_related[200:250,18:23]
tail(potomac_for_related[,10:23])


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#   #All Potomac -  All Within-Site Reps    ###########
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

## Keep only sample label and the genotype data
## Get rid of the MLG designations and the headers

## To get rid of headers, you first convert the dataframe to a matrix, then strip out the header, then convert the matrix back to a dataframe in which the variables are not factors

allpotomac.related.dataframe <-potomac_for_related %>% 
  dplyr::select(-c(IDName,NewPop,TideCode, Clone.ID.2018))%>% 
  as.matrix(.) %>% 
  matrix(., ncol = ncol(.), dimnames = NULL) %>% #remove variable names
  as.data.frame(., stringsAsFactors=FALSE) #convert back to a data.frame

names(allpotomac.related.dataframe)
head(allpotomac.related.dataframe)

## Write out data ----
## Write out tab delimited version of the data in a form that related can use so you don't have to recreate it again.

write.table(allpotomac.related.dataframe, here("processed.data","related", "allsamples","relateddata.allsamples.potomac.txt"))

###---Convert data into format used by related ----
allpotomac.relateddata <- readgenotypedata(allpotomac.related.dataframe)

###---Compare relatedness estimators ----
#compareestimators(allpotomac.relateddata, 100)

###---Estimate relatedness ---####
coancestry.results.allsamples <- coancestry(allpotomac.relateddata$gdata, wang=1)

## you could calculate the other estimators as well if you want but this code is only doing wang's r

####---Save output from wang estimator to Rdata ----

save(coancestry.results.allsamples,file = here("processed.data","related", "allsamples","coancestry.results.allsamples.RDATA"))

#XXXXXXXXXXXXXXXXXXX
#No Reps   #########
#XXXXXXXXXXXXXXXXXXX

##keep only one rep of each MLG at each site and after filtering, Keep only the concatenated sample label and the genotype data.

#As with the full data set, get rid of headers by converting to a matrix, removing headers and converting back to a data.frame.

noreps.related.dataframe <-potomac_for_related %>% 
  dplyr::distinct(NewPop, Clone.ID.2018, .keep_all = TRUE) %>% 
  dplyr::select(-c(IDName,NewPop,TideCode, Clone.ID.2018)) %>% 
  as.matrix(.) %>% 
  matrix(., ncol = ncol(.), dimnames = NULL) %>% #remove variable names
  as.data.frame(., stringsAsFactors=FALSE) #convert back to a data.frame

head(noreps.related.dataframe)

## Write out data ----

##Save the data set that related will use so you do not have to recreate it over and over if you need to redo the analysis

write.table(noreps.related.dataframe, here("processed.data","related", "noreps","relateddata.noreps.potomac..txt"))

####---Convert data into format used by related----
noreps_potomac.relateddata <- readgenotypedata(noreps.related.dataframe)

###---Compare relatedness estimators----
#compareestimators(noreps_potomac.relateddata, 100) # ran this to check which estimator to use.

###---Estimate relatedness using wangs  estimator----
coancestry.results.noreps <- coancestry(noreps_potomac.relateddata$gdata, wang=1)

## you could calculate the other estimators as well if you want but this code is only doing wang's r

###---Write output to a dataframe ---- 

#Save the relatedness dataframe
save(coancestry.results.noreps,file = here("processed.data","related", "noreps","coancestry.results.noreps.RDATA"))

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# grouprel analysis of significance of relatedness in groups ####
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#---Test whether samples or MLGs within tidal and nontidal environments and in the whole river are more closely related than expected ---# 

#These lines are commented out because the runs take a very long time - ~11-12 hours for the full data set.  You need to make sure you really want to run them - and if you do, just uncomment the lines.

#The results of grouprel are written to files with the generic names in the working directory - observed-r.csv, expectedrel.csv, and Rplot.pdf. So you will want to either change working directories before you run the second one, or rename the output files so they do not get overwritten.

#AB=Tidal, AA=Nontidal
# testgroup_mlgs_related<-related::grouprel(genotypes = noreps_potomac.relateddata$gdata, 
#                                      estimatorname = "wang",
#                                      usedgroups = "all", iterations = 1000)
# 
# 
# testgroup_allsamples.related<-related::grouprel(genotypes = allpotomac.relateddata$gdata,
#                                                 estimatorname = "wang",
#                                                 usedgroups = "all", iterations = 1000)
# 

