## Author: Maile Neel
## 2024

## This script uses functions from adegnet, pegas, hierfstat, and mmod to summarize data in genind objects created in the 01_Create_adegenet_geninds.R script that were saved to the file list.of.genind.objects.rds in the processed.data folder. If that file does not exist, you need to create it before proceeding here. If it does exist proceed to read it in.

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Install and Load Packages & Source Functions ----
#XXXXXXXXXXXXXXXXXXXXxxxXXXXXXXXXXXXXXXXXXXXXXXXXXXX

## This script uses a version of mmod that removes missing data from GST bootstrap calculations. Install it first before running pacman

# devtools::install_github("dwinter/mmod", ref="bs_fix")

if (!require("pacman")) install.packages("pacman")
pacman::p_load(here,
               tidyverse,
               devtools,
               adegenet,
               pegas,
               mmod,
               hierfstat,
               reshape2)

if(!exists("calculate.popgen.diversity.withinsite", mode="function")) 
  source(here::here("functions", "calculate.popgen.diversity.withinsite.r"))

if(!exists("calculate.popgen.diversity.amongsite", mode="function")) 
  source(here::here("functions", "calculate.popgen.diversity.amongsite.r"))

#XXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Read in data objects ----
#XXXXXXXXXXXXXXXXXXXXXXXXXXXX

## Read in the rds file that contains the genind objects created previously, and filter for only the geninds based on the unique MLGs.

## If the rds file is not in the processed.data folder, create it with the 01.PR_Manuscript_adegenet_Create_Dataframes.r script.

list.of.genind.objects.for.diversity <- 
  readRDS(here("processed.data",
               "list.of.genind.objects.rds")) %>% 
  imap(~if (grepl("noreps", .y)) .x else NULL) %>%
  compact() 

names(list.of.genind.objects.for.diversity) <-  sub("genind.", "", names(list.of.genind.objects.for.diversity))


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#  Within-Site Calculations ----
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#Use custom function to summarize the within site statistics for the overall data set. The statistics are calculated using adegenet and hierfstat. See the custom function for details on which statistic uses which package.

## min.alleles is needed by hierfstat to calculate the number of alleles rarefied to the smallest sample at any site. We include that as an argument in our function
  
within.site.diversity<- calculate.popgen.diversity.withinsite(list.of.genind.objects.for.diversity$noreps.potomac) 
  

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#  Among-Site and Tidal Regime Calculations ----
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#++++++++++++++++++++++++++
## Global Calculations ----
#++++++++++++++++++++++++++

## Calculate the statistics for four data sets included in the list of genind objects.

## use mmod package to get Hedrick’s G”st and Jost’s D
list.of.global.diff.stats <- list.of.genind.objects.for.diversity %>% 
  map(., ~mmod::diff_stats(.x))


combined_lists.global.diff.stats <- list()

# Loop through each named list
for (name in names(list.of.global.diff.stats)) {
  # Extract the global dataframe
  global_df <- as.data.frame(t(list.of.global.diff.stats[[name]]$global))
  
  # Add a new column with the higher-level list name
  global_df$list_name <- name
  
  # Append the dataframe to the combined list
  combined_lists.global.diff.stats <- append(combined_lists.global.diff.stats, list(global_df))
}

# Combine all dataframes into one
combined.global.diff.stats.df <- bind_rows(combined_lists.global.diff.stats) %>% 
  select(list_name,Gprime_st,D_het)

## use hierfstat to get Weir and Cockeram's FST

list.of.global.wc.fstats <- list.of.genind.objects.for.diversity %>% 
  purrr::map(., ~genind2hierfstat(.x)) %>% 
  map(., ~wc(.x))

combined_lists.wc.stats <- list()

# Loop through each named list
for (name in names(list.of.global.wc.fstats)) {
  # Extract the global dataframe
  global_df <- list.of.global.wc.fstats[[name]]$FST 
    names(global_df)<-"FST"
  
  # Add a new column with the higher-level list name
  global_df$list_name <- name
  
  # Append the dataframe to the combined list
  combined_lists.wc.stats <- append(combined_lists.wc.stats, list(global_df))
}

# Combine all wc estimate dataframes into one

combined_lists.wc.stats.df <- bind_rows(combined_lists.wc.stats)

## Join the FST values with mmod values by dataset.

final.all.global.among.population.measures <- combined_lists.wc.stats.df %>% 
  left_join(combined.global.diff.stats.df, by = "list_name") %>% 
  dplyr::relocate(list_name, .before=FST) 

#+++++++++++++++++++++++++++++++++++++++++++++++++++
### Write global among-site statistics to csv file ----
#+++++++++++++++++++++++++++++++++++++++++++++++++++

write_csv(final.all.global.among.population.measures, 
          here("processed.data", 
               "final.all.global.among.population.measures.csv"))


### Bootstrap Confidence Intervals ----

wc_bootstrap <- function(genind) {
  wc.boot<-boot.vc(genind2hierfstat(genind)[,1],
                genind2hierfstat(genind)[,-1],
                nboot=1000, 
                quant=c(0.025,0.5,0.975))
  wc.confidence<-wc.boot$ci[,2]
}

FST_confidence.intervals<-map(list.of.genind.objects.for.diversity, wc_bootstrap) %>% 
as.data.frame() %>% 
  rownames_to_column()
 
## Function for calculating bootstraps with MMOD

mmod_bootstrap <- function(genind, bs) {
  bs <- chao_bootstrap(genind, nreps=1000) #generate bootstraps
 
  obs.D <- D_Jost(genind) 
  bs_D <- summarise_bootstrap(bs, D_Jost) 
  bias_D <- bs_D$summary.global.het[1]- obs.D$global.het 
  D_ci <- bs_D$summary.global.het- bias_D
  
  obs.Gst <- Gst_Hedrick(genind) 
  bs_Gst <- summarise_bootstrap(bs, Gst_Hedrick) 
  bias_Gst <- bs_Gst$summary.global[1]- obs.Gst$global 
  Gst_CI <-bs_Gst$summary.global.het- bias_Gst
  
mmod_CIs<-as.data.frame(rbind(bs_D$summary.global.het-bias_D,bs_Gst$summary.global.het- bias_Gst)) %>% 
    mutate(Metric = c("Dest", "GST")) %>% 
    relocate(Metric)
  
  return(mmod_CIs)
}

## map function onto the list of geninds.
MMOD_confidence.intervals<-map(list.of.genind.objects.for.diversity, mmod_bootstrap, bs)

MMOD_confidence.intervals.df <- map_dfr(MMOD_confidence.intervals, 
                                        ~mutate(.x, source = cur_group_id()), 
                                        .id = "source") %>% 
  relocate(source)

# Save the two confidence interval objects for later use


#### Write confidence interval data to csv file ----


write_csv(MMOD_confidence.intervals.df, 
          here("processed.data","MMOD_confidence.intervals.df.csv"))

write_csv(FST_confidence.intervals, 
          here("processed.data","FST_confidence.intervals.csv"))

#++++++++++++++++++++++++++++++++++
## Pairwise Among-Site Calculations ----
#++++++++++++++++++++++++++++++++++

list.of.pairwise.wc.fstats <- list.of.genind.objects.for.diversity %>% 
  map(., ~hierfstat::pairwise.WCfst(.x)) %>% 
  map(., ~calculate.popgen.diversity.amongsite(.x, metric = "FST"))

list.of.pairwise.GST.Hedrick <- list.of.genind.objects.for.diversity %>% 
  map(., ~mmod::pairwise_Gst_Hedrick(.x)) %>% 
  map(., ~calculate.popgen.diversity.amongsite(.x, metric = "Gst.Hed"))

list.of.pairwise.JostsD <- list.of.genind.objects.for.diversity %>% 
  map(., ~mmod::pairwise_D(.x)) %>% 
  map(., ~calculate.popgen.diversity.amongsite(.x, metric = "Josts.D"))


## Join the dataframes for each among-site metric into one

## Function to do the joining
join_dataframes <- function(df_list) {
  Reduce(function(x, y) full_join(x, y, by = "Pop"), df_list)
}

##Loop the function over our lists

combined.among.site.lists <- list()

for (name in names(list.of.pairwise.wc.fstats)) {
  combined.among.site.lists[[name]] <- 
    join_dataframes(list(list.of.pairwise.wc.fstats[[name]], 
                         list.of.pairwise.GST.Hedrick[[name]], 
                         list.of.pairwise.JostsD[[name]]))
  }

## Create one data frame for the among-site statistics calculated from sites only within each tidal regime

stats.within.regime <- rbind(combined.among.site.lists$noreps.potomac.NonTidal,combined.among.site.lists$noreps.potomac.Tidal) %>% 
  dplyr::rename(FST.within = FST, Gst.Hed.within = Gst.Hed, Josts.D.within = Josts.D)

## Join the within regime values with the values from the whole river and with the within population statistics into one dataframe

all.popgen.stats <- within.site.diversity %>% 
  dplyr::left_join(combined.among.site.lists$noreps.potomac, by = "Pop") %>% 
  dplyr::left_join(stats.within.regime, by = "Pop")
  

## To return the OrderPop and Tide variables to the corresponding site data, we need to pull those values out of the @other slot of one genind and get the unique values and then join the values to each of the dataframes in the final list.

site.labels <-  as.data.frame(list.of.genind.objects.for.diversity$noreps.potomac@strata) %>% 
  dplyr::group_by(NewPop) %>% 
  dplyr::summarize(OrderPop = first(OrderPop),
                   Tide = first(Tide)) 


final.site.level.popgen.stats.with.labels <- site.labels %>% 
  dplyr::left_join(all.popgen.stats, by = c("NewPop" = "Pop"))

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Write site level statistics to csv file ----
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

write_csv(final.site.level.popgen.stats.with.labels, 
          here("processed.data", "final.site.level.popgen.stats.with.labels.csv")) 
