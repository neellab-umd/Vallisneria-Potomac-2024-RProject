## Author: Maile Neel

## Usage ----

## The script is designed to be run within the R Project Calm After the Storm.Rproj.  It uses cost distances among sample sites in the Hudson River that are read in from the csv file CATS_Site_Cost_Distances_combined.csv in the output.data directory if that script has been run.

## The distances were calculated in the script D3A_Calculate_Site_Water_Distances.r using a transition matrix created in the script CATS_D2_gdistance_surface_code_6_2021.r  


## Load packages ----
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,
               here,
               ggthemes,
               update = FALSE)


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## 1. Read in distance data ----  
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## In this case the points are the centroids of the sample locations used in the CATS manuscript. To make summarizing and plotting easier, replace the year 2014 from CPN with 2015.

site.rescost.distlist_combined <- read_csv(here::here("processed.data",
                                                      "Site_Cost_Distances.csv")) 
  
#+++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++
## 2. Summarize and Plot Distances Among Sites 
#+++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++

## Summarize distances among sites by the "From" sample. Summarizing on the To Sample column would give you the exact same result because each sample is present in each column.  So here, for each "From" sample, we get the distances to all other samples, but we do not retain which "To" sample is which

## Plot intersite distances from each pop to all other pops within tidal regime.

ggplot.site.rescost.distlist_combined <- 
  site.rescost.distlist_combined %>% 
  dplyr::mutate(TideClass = case_when(FromTide == ToTide ~"Same",
                                    FromTide !=ToTide ~"Different"),
                FromTide = factor(FromTide, levels = c("Nontidal", "Tidal")),
                FromPop = fct_reorder(FromPop,FromOrderPop)) %>% 
  dplyr::distinct(FromPop,ToPop, .keep_all = TRUE) %>%
  dplyr::filter(TideClass == "Same") %>% 
  dplyr::group_by(FromPop) %>% 
  ggplot(aes(x = FromPop, y = WaterDist))+
  stat_summary(fun = "mean", shape=17, color = "gray80",size=.5)+
  geom_jitter(width = 0.1, height = 0, size = 1.5)+
  labs(y = "Distances to Other Sites (km)", x = "Sample Site")+
  facet_grid(~FromTide, scales="free_x", space = "free_x")+
  theme_stata(scheme = "s1mono")+
  theme(axis.title = element_text(size = 13),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .25),
        axis.text.y = element_text(hjust = .5))


ggsave(here("figures", "Appendix2.WaterDistances.Among.Sites.By.Tide.tiff"),
       plot = last_plot(),
       width = 16, height = 10, units = c("cm"),
       dpi = 350,
       limitsize = TRUE)


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## 2. Summarize Distances Among Sites 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## Summarize the minimum, median, and maximum distances for site using the FromPop variable. 

Site.Water.Dist.Summary.Within.tide<-site.rescost.distlist_combined %>%
  dplyr::mutate(TideClass = case_when(FromTide == ToTide ~"Same",
                                      FromTide !=ToTide ~"Different"),
                FromTide = factor(FromTide, levels = c("Nontidal", "Tidal"))) %>% 
  dplyr::group_by(FromPop,FromOrderPop) %>%
  dplyr::filter(TideClass == "Same")%>%
  dplyr::summarise(Tide = unique(FromTide), MedDistance.within.tide=median(WaterDist), MeanDistance.within.tide=mean(WaterDist),MinDistance.within.tide=min(WaterDist),MaxDistance.within.tide = max(WaterDist))

Min.Water.Distances.by.tide<-Site.Water.Dist.Summary.Within.tide %>% 
  group_by(Tide) %>% 
  summarize_at(vars(MinDistance.within.tide), list(mean=mean, med=median,min=min,max=max))

Mean.Water.Distances.by.tide<-Site.Water.Dist.Summary.Within.tide %>% 
  group_by(Tide) %>% 
  summarize_at(vars(MeanDistance.within.tide), list(mean=mean, med=median,min=min,max=max))

Max.Water.Distances.by.tide<-Site.Water.Dist.Summary.Within.tide %>% 
  group_by(Tide) %>% 
  summarize_at(vars(MaxDistance.within.tide), list(mean=mean, med=median,min=min,max=max))

Within.Tide.Water.Dist.Summary <- site.rescost.distlist_combined %>%
  dplyr::mutate(TideClass = case_when(FromTide == ToTide ~"Same",
                                      FromTide !=ToTide ~"Different"),
                FromTide = factor(FromTide, levels = c("Nontidal", "Tidal"))) %>% 
  dplyr::filter(TideClass == "Same")%>%
  dplyr::group_by(FromTide) %>%
  dplyr::summarise(Tide = unique(FromTide), MedDistance.within.tide=median(WaterDist), MeanDistance.within.tide=mean(WaterDist),MinDistance.within.tide=min(WaterDist),MaxDistance.within.tide = max(WaterDist))

Tide.Min.Water.Distances.by.tide<-Site.Water.Dist.Summary.Within.tide %>% 
  group_by(Tide) %>% 
  summarize_at(vars(MinDistance.within.tide), list(mean=mean, med=median,min=min,max=max))

Tide.Mean.Water.Distances.by.tide<-Site.Water.Dist.Summary.Within.tide %>% 
  group_by(Tide) %>% 
  summarize_at(vars(MeanDistance.within.tide), list(mean=mean, med=median,min=min,max=max))

Tide.Max.Water.Distances.by.tide<-Site.Water.Dist.Summary.Within.tide %>% 
  group_by(Tide) %>% 
  summarize_at(vars(MaxDistance.within.tide), list(mean=mean, med=median,min=min,max=max))
