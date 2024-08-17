## Author: Maile Neel
## 2024
 
## Usage --

## This script generates summaries of distances among individuals for the Potomac River that were calculated in 08A_Calculate_Sample_Distances.R.  That script creates a csv file called allpotomac.euclidean.and.cost.distances.among.samples.csv that is in the processed.data folder.  The csv can be read in here if it is not already loaded.  

## If allpotomac.euclidean.and.cost.distances.among.samples.csv is not in the processed.data directory, you need to run 08A_Calculate_Sample_Distances.R.

## The nearest neighbor and maximum distances among samples within sites calculated here are presented in an appendix in the manuscript.

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Install and Load Packages  ----
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,
               magrittr,
               here,
               ggpubr,
               scales,
               ggforce,
               ggthemes,
               update = FALSE
               )

#XXXXXXXXXXXXXXXXXXXXXXX
# Read in csv File  ----
#XXXXXXXXXXXXXXXXXXXXXXX

## Read in csv of pairwise distances among samples created in the 08A_Calculate_Sample_Distances.R.if the object Ind.Euclidean.and.Cost.Dist.List is not in your environment

Ind.Euclidean.and.Cost.Dist.List  <- 
  read_csv(here("processed.data", 
                "allpotomac.euclidean.and.cost.distances.among.samples.csv"))


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Summarize Distances Among Samples with Sites ----
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

### Summarize to Sample ----

## Summarize the minimum, median, and maximum distances for each sample (using the FromSample variable) in each site (using the FromPop variable). For each sample, we use the distances to all other samples, but we do not retain which "to" sample is which. Summarizing on the To Sample column would give you the exact same result because each sample is present in each column.

## FromMLG and FromTide are included in the group_by because they are needed for further summaries. They do not select any further than FromPop and FromSample in this query.

distsumall_by.sample.within.pop  <- Ind.Euclidean.and.Cost.Dist.List %>%
  dplyr::filter(SameSite ==  "Same")%>%
  dplyr::group_by(FromSample, FromPop, FromMLG, FromTide, FromOrderPop) %>%
  dplyr::summarise(MaxDistance = max(FinalDist), 
                   MinDistance = min(FinalDist),
                   MeanDistance = mean(FinalDist),
                   MedDistance = median(FinalDist))

## You can write this intermediate summary in which each sample is one line out to a csv if you want to.

#write_csv(distsumall_by.sample.within.pop, file = here("processed.data", "allpotomac.euclidean.and.cost.distsumall.by.Sample.within.Pop.csv"))

### Summarize Samples to Site ----

## Further summarize the nearest neighbor distance and maximum distance among samples within sites.  For each of those, we calculate minimum, median, and maximum distances among all samples in each site

distsumall_by.pop  <- distsumall_by.sample.within.pop %>%
  dplyr::group_by(FromPop, FromTide, FromOrderPop) %>%
  dplyr::summarise(n = n(),
                   MinMinDistance = min(MinDistance),
                   MeanMinDistance = 
                     mean(MinDistance),
                   MedMinDistance = median(MinDistance),
                   MaxMinDistance = max(MinDistance), 
                   MeanMeanDistance = mean(MeanDistance),
                   MinMaxDistance = min(MaxDistance),
                   MeanMaxDistance = 
                     mean(MaxDistance),
                   MedMaxDistance = median(MaxDistance),
                   MaxMaxDistance = max(MaxDistance)) 

distsumall_by.pop %>% 
mutate(across(c(MinMinDistance:MaxMaxDistance), \(x) round(x, 1))) %>% 
  arrange(FromOrderPop) %>% 
  write_csv(., file = here("processed.data",
                           "allpotomac.euclidean.and.cost.distsumall.allsamples.within.sites.csv"))

### Summarize Sites to Tidal Regime ----

## This code creates the summary of site level values of distances among all samples within sites for all sites in each tidal regime.

distsumall_by.tide <- distsumall_by.pop %>%
  dplyr::group_by(FromTide) %>%
  dplyr::summarise(n = n(),
                   MinofMinDistance = min(MeanMinDistance),
                   MeanofMinDistance = mean(MeanMinDistance),
                   MaxofMinDistance = max(MeanMinDistance), 
                   MedianofMeanDistance = median(MeanMeanDistance),
                   MeanofMeanDistance = mean(MeanMeanDistance),
                   MinofMaxDistance = min(MeanMaxDistance),
                   MeanofMaxDistance = mean(MeanMaxDistance),
                   MaxofMaxDistance = max(MeanMaxDistance))

  
## Plot distances among samples within sites.
## First prep the data set for plotting by filtering for only samples from the same site and creating factors to enforce plotting order.  Use dplyr::reframe() instead of summarise() because summarise() can't be used when you have more than one row per grouped condition and there are multiple cases of distances between the same FromMLG and ToMLG. 

ggplot.Within.Pop.Sample.Distances <-
  Ind.Euclidean.and.Cost.Dist.List %>% 
  dplyr::filter(FromPop == ToPop) %>% 
  dplyr::mutate(FromTide = 
                  factor(FromTide, levels = c("NonTidal", "Tidal")),
                FromOrderPop = factor(FromOrderPop)
                ) %>% 
  dplyr::distinct(FromSample,ToSample, .keep_all = TRUE) %>% 
  dplyr::group_by(FromOrderPop,FromMLG,ToMLG,FromTide) %>%
  dplyr::reframe(FinalDist = identity(FinalDist)) %>% 

## Make the ggplot
  ggplot(aes(FromOrderPop, FinalDist))+
  geom_violin(fill = NA, linewidth = .85) +
  ggforce::geom_sina(size = .5) +
  labs(y = "Distances Among Samples (m)", x = "Sample Site") +
  facet_grid(~FromTide, scales = "free_x", space = "free_x") +
  geom_point(data = distsumall_by.pop, aes(x = as.factor(FromOrderPop), y=MinMinDistance), 
             color = "gray80", fill = "gray80",
             shape = 25, size = 1.5) +
  geom_point(data = distsumall_by.pop, aes(x = as.factor(FromOrderPop), y=MeanMeanDistance),
             color = "gray80",fill = "gray80",
             shape = 24, size = 1.5) +
  ggthemes::theme_stata(scheme = "s1mono") +
  theme(axis.title = element_text(size = 13),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        axis.text.y = element_text(hjust = .5))

ggsave(here("figures", "Fig.Appendix.2.Distances.Among.Samples.By.Pop.And.Tide.tiff"),
            plot = ggplot.Within.Pop.Sample.Distances,
            width = 15, height = 10, units = c("cm"),
            dpi = 350,
            limitsize = TRUE)

#Summarize Distances of MultiSample Clones ----

CloneDist <- Ind.Euclidean.and.Cost.Dist.List %>% 
  dplyr::filter(ToMLG == FromMLG & ToSample != FromSample) %>% 
  dplyr::group_by(ToMLG, ToTide) %>% 
  dplyr::summarize(mean = mean(FinalDist)/1000,
            max = max(FinalDist)/1000,
            num_sites = n_distinct(ToPop)) %>% 
  arrange(ToTide, desc(max))

write_csv(CloneDist, file = here("processed.data",
                                 "multi.site.mlg.distances.csv"))
