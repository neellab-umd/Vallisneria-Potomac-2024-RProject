## Author: Maile Neel
## 2024

## Usage --

## This script is set up to use within the RProject Potomac.Manuscript.2024.R.Project.Rproj. It tests IBD using Mantel Test of correlations between distances among sites that were calculated in 09A_Calculate_Site_Water_Distances.r and genetic distance calculated here using one of the genpop object created in the script 01_Create_adegenet_geninds.R that was written to ./processed.date/list.of.genpop.objects.rds


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Install and Load Packages  ----
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

if (!require(pacman)) install.packages("pacman")
pacman::p_load(adegenet,
               ade4,
               tidyverse,
               here,
               magrittr,
               update = FALSE)


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Read in Required Files  ----
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXX

## The code below requires the gen.pop objects created in the adegenet 2018_PR_Manuscript_adegenet_create_dataframes.R script.  If the object below is not in the processed.data folder, run script first to generate the object.

list.of.genpop.objects <- readRDS(here("processed.data","list.of.genpop.objects.rds"))

## It also requires the matrices with distances among sites created using the scripts 09a_Calculated_Site_Water_Distances.R.  If the rds files containing these matrices are not in the processed.data folder, that script needs to be run first.

Site.Final.dist <- readRDS(here("processed.data", "Full_PR_Site_Cost_and_Euclidean_Distances_combined_Matrix.rds"))

Site.Final.dist.tidal <- readRDS(here("processed.data", "Tidal_PR_Site_Cost_and_Euclidean_Distances_combined_Matrix.rds"))

Site.Final.dist.nontidal <- readRDS(here("processed.data", "Nontidal_PR_Site_Cost_and_Euclidean_Distances_combined_Matrix.rds"))

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
## IBD Tests For All Samples ####
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

### All Sites ----
#### Calculate genetic distance ----

## The genetic distance options in adegenet are as follows:
## 1) Nei's; 2) Edward's chord distance; #3) Coancestralilty; 4) classic Euclidean; 5) Absolute (Provesti's)

## Using the gen.pop objects in the list read in above, I used method 2 -Angular distance or Edwards' Chord distance - because it is a Euclidean distance measure.

Distgen_all <- adegenet::dist.genpop(list.of.genpop.objects$genpop.allpotomac,
                           method = 2, upper = TRUE)

#### Test IBD ----
 
## Using the distance matrix from the gdistance cost analysis
## Need to make sure the output matrix from gdistance is a dist object (using as.dist)

ibd_allsamp_allsites_cost <- ade4::mantel.randtest(Distgen_all,
                                                   as.dist(Site.Final.dist))

ibd_allsamp_allsites_cost
plot(ibd_allsamp_allsites_cost)

####Basic IBD plot ----
plot(as.dist(Site.Final.dist/1000), Distgen_all, 
     xlab = "Cost Distance (km)", 
     ylab = "Edward's Chord Genetic Distance")
abline(lm(as.vector(Distgen_all) ~as.vector(Site.Final.dist/1000)), 
       col = "black",lty = 1, lwd = 2)

### NonTidal Sites ----

#### Calculate genetic distance ----

Distgen_all_NonTidal <- adegenet::dist.genpop(list.of.genpop.objects$genpop.allpotomac.NonTidal,
                                    method=2)

#### Test IBD ----

ibd_allsamp_NonTidal <- ade4::mantel.randtest(Distgen_all_NonTidal,
                                              as.dist(Site.Final.dist.nontidal))

ibd_allsamp_NonTidal
plot(ibd_allsamp_NonTidal)

####Basic IBD plot ----

plot(as.dist(Site.Final.dist.nontidal/1000), Distgen_all_NonTidal, 
     xlab = "Cost Distance (m)", 
     ylab = "Edward's Chord Distance")
abline(lm(as.vector(Distgen_all_NonTidal) ~as.vector(Site.Final.dist.nontidal/1000)), 
       col = "black",lty = 1, lwd = 2)

###Tidal Sites ----

#### Calculate Genetic distance ----
Distgen_all_Tidal <- adegenet::dist.genpop(list.of.genpop.objects$genpop.allpotomac.Tidal,
                                 method = 2)

#### Test IBD ----

ibd_allsamp_Tidal <- ade4::mantel.randtest(Distgen_all_Tidal,
                                           as.dist(Site.Final.dist.tidal))

ibd_allsamp_Tidal
plot(ibd_allsamp_Tidal)

#### Basic IBD plot ----
plot(as.dist(Site.Final.dist.tidal/1000), Distgen_all_Tidal, 
     xlab = "Cost Distance (m)", 
     ylab = "Edward's Chord Distance")
abline(lm(as.vector(Distgen_all_Tidal) ~as.vector(Site.Final.dist.tidal/1000)), 
       col = "black",lty = 1, lwd = 2)


#+++++++++++++++++++++++++++++++++++++++++++++++++
## IBD Tests For Each MLG In Each Population####
#+++++++++++++++++++++++++++++++++++++++++++++++++

###All Sites ----

#### Calculate Genetic distance ----

Distgen_noreps <- adegenet::dist.genpop(list.of.genpop.objects$genpop.noreps.potomac,
                                        method=2)

#### Test IBD ----

#This is using the distance matrix from the gdistance cost analysis in which the distance is the greater of either cost or Euclidean.

ibd_noreps_allsites_cost <- ade4::mantel.randtest(Distgen_noreps,
                                                  as.dist(Site.Final.dist))

ibd_noreps_allsites_cost
plot(ibd_noreps_allsites_cost)

####Basic IBD plot ----

plot(as.dist(Site.Final.dist), Distgen_noreps)
abline(lm(as.vector(Distgen_noreps)~as.vector(Site.Final.dist)), 
       col="red",lty=2)


### NonTidal Sites only ----

#### Calculate Genetic distance ----

Distgen_noreps_NonTidal <- 
  adegenet::dist.genpop(list.of.genpop.objects$genpop.noreps.potomac.NonTidal,
                                     method = 2)

#### Test IBD ----

ibd_noreps_NonTidal <- ade4::mantel.randtest(Distgen_noreps_NonTidal,
                                             as.dist(Site.Final.dist.nontidal))

plot(ibd_noreps_NonTidal)

#### Basic IBD plot ----
plot(as.dist(Site.Final.dist.nontidal), Distgen_noreps_NonTidal)
abline(lm(as.vector(Distgen_noreps_NonTidal) ~as.vector(Site.Final.dist.nontidal)), 
       col = "red",lty = 2)

### Tidal Sites only ----

#### Calculate Genetic distance ----

Distgen_noreps_Tidal <- adegenet::dist.genpop(list.of.genpop.objects$genpop.noreps.potomac.Tidal,
                                  method=2)

#### Test IBD ----

ibd_noreps_Tidal <- ade4::mantel.randtest(Distgen_noreps_Tidal,
                                          as.dist(Site.Final.dist.tidal),nrepet=999)
ibd_noreps_Tidal


#### Basic IBD plot ----

plot(as.dist(Site.Final.dist.tidal), Distgen_noreps_Tidal)
abline(lm(as.vector(Distgen_noreps_Tidal) ~as.vector(Site.Final.dist.tidal)), 
       col = "red",lty = 2)

###Gather up all the IBD tests into a list and pull out the 
###observed correlation  and p values.
Mantel.results <- mget(ls(pattern = "ibd_"))

Mantel.results.summary<-data.frame(Mantel.results.summary = 
                        sapply(Mantel.results, function(x) x$obs), 
                        p.value = sapply(Mantel.results, function(x) x$pvalue))

write.csv(Mantel.results.summary, file = "./processed.data/Potomac.Mantel.Test.Summary.csv")

## Data wrangling for ggplot ----
## Prepare for creating a ggplot of genetic by geographic distance by pulling both distances into one dataframe, mutating from distances to numeric, and pivoting to put the all reps distances and unique mlg distances into one column so they can be plotted in ggplot as two groups in the same variable

### Full River ----

IBD_full_river_df <- tibble(Distgen_all,Distgen_noreps,
                             Site.Final.dist) %>% 
 mutate_all(as.numeric) %>% 
 pivot_longer(cols = c(Distgen_all,Distgen_noreps),
               names_to = "Distance_Type", 
               values_to = "Genetic_Distance") %>% 
  mutate(Set = "Full River") %>% 
  rename(Geographic_Distance = Site.Final.dist)

### Nontidal ----

IBD_nontidal_df <- tibble(Distgen_all_NonTidal,Distgen_noreps_NonTidal,
                                                 Site.Final.dist.nontidal)%>% 
  mutate_all(as.numeric) %>% 
  pivot_longer(cols = c(Distgen_all_NonTidal,Distgen_noreps_NonTidal),
               names_to = "Distance_Type", 
               values_to = "Genetic_Distance") %>% 
  mutate(Distance_Type = str_replace(Distance_Type, "_NonTidal", ""),
         Set = "Nontidal") %>% 
  rename(Geographic_Distance = Site.Final.dist.nontidal)

### Tidal ----

IBD_tidal_df <- tibble(Distgen_all_Tidal,Distgen_noreps_Tidal,
                                              Site.Final.dist.tidal) %>% 
  mutate_all(as.numeric) %>% 
  pivot_longer(cols = c(Distgen_all_Tidal,Distgen_noreps_Tidal),
               names_to = "Distance_Type", 
               values_to = "Genetic_Distance") %>% 
  mutate(Distance_Type = str_replace(Distance_Type, "_Tidal", ""),
         Set = "Tidal") %>% 
  rename(Geographic_Distance = Site.Final.dist.tidal)

### Combine the full river, nontidal and tidal data ----

IBD_All_Sets <- rbind(IBD_full_river_df,IBD_nontidal_df,IBD_tidal_df)

## Create ggplots of ibd ----

theme_set(theme_classic())

### Multiplot by River Reach ----

## Create individual plots for 1) the full river, 2) tidal, and 3) nontidal coded by whether the genetic distance is based on all samples or just unique MLGs. 

#### Full River ----
  IBD_Plot_FullRiver <- IBD_All_Sets %>% 
    filter(Set == "Full River" ) %>% 
    ggplot() +
  geom_point(aes(x = Geographic_Distance/1000, 
                 y = Genetic_Distance, color = Distance_Type)) +
  geom_smooth(aes(x = Geographic_Distance/1000, 
                  y = Genetic_Distance, color = Distance_Type),
              method = "lm", se = FALSE) +
  annotate("text", label = "Full River", x = 0, y = 0.62, hjust = 0) +
  ylab(label = "Edward's Chord Genetic Distance") +
  xlab(label = "Distance Through Water") +
  xlim(c(0, max(IBD_All_Sets$Geographic_Distance)/1000)) +
  ylim(c(0, max(IBD_All_Sets$Genetic_Distance))) +
  theme(legend.position = "none") +
  scale_color_grey()

#### Tidal ----

  IBD_Plot_Tidal <- IBD_All_Sets %>% 
    filter(Set == "Tidal" ) %>% 
    ggplot() +
    geom_point(aes(x = Geographic_Distance/1000, 
                   y = Genetic_Distance, color = Distance_Type)) +
    geom_smooth(aes(x = Geographic_Distance/1000, 
                    y = Genetic_Distance, color = Distance_Type),
                method = "lm", se = FALSE) +
    annotate("text", label = "Tidal", x = 0, y = 0.62, hjust = 0) +
    ylab(label = "Edward's Chord Genetic Distance") +
    xlab(label = "Distance Through Water") +
    xlim(c(0, max(IBD_All_Sets$Geographic_Distance)/1000)) +
    ylim(c(0, max(IBD_All_Sets$Genetic_Distance))) +
    theme(legend.position = "inside", 
          legend.position.inside = c(0.8, 0.4),
          legend.title = element_blank()) +
    scale_color_grey(labels = c("Distgen_all_Tidal" = "All Samples", 
                                "Distgen_noreps_Tidal" = "Unique MLGs"))

#### Nontidal ----

  IBD_Plot_Nontidal <- IBD_All_Sets %>% 
    filter(Set == "Nontidal" ) %>% 
    ggplot() +
    geom_point(aes(x = Geographic_Distance/1000, 
                   y = Genetic_Distance, color = Distance_Type)) +
    geom_smooth(aes(x = Geographic_Distance/1000, 
                    y = Genetic_Distance, color = Distance_Type),
                method = "lm", se = FALSE) +
    annotate("text", label = "Nontidal", x = 0, y = 0.62, hjust = 0) +
    ylab(label = "Edward's Chord Genetic Distance") +
    xlab(label = "Distance Through Water") +
    xlim(c(0, max(IBD_All_Sets$Geographic_Distance)/1000)) +
    ylim(c(0, max(IBD_All_Sets$Genetic_Distance))) +
    theme(legend.position = "none") +
    scale_color_grey()
  
#### Combined Plot ----

IBD_Plot_Combined_by_Tide <- IBD_Plot_FullRiver/IBD_Plot_Nontidal/IBD_Plot_Tidal+
  patchwork::plot_layout(axis_titles = "collect")

ggsave(plot = IBD_Plot_Combined_by_Tide, here('figures','IBD_Plot_Combined_by_Tide.png'), width = 10, height = 15, units = "cm", dpi=400)

### Multiplot by Distance Type

## Create individual plots for 1) genetic distance calculated by using all samples and 2) genetic distance based only on unique MLGs 

#### All Samples ----

IBD_Plot_All_Samples <- IBD_All_Sets %>% 
  filter(Distance_Type == "Distgen_all") %>% 
  ggplot() +
  geom_point(aes(x = Geographic_Distance/1000, 
                 y = Genetic_Distance, color = Set)) +
  geom_smooth(aes(x = Geographic_Distance/1000, 
                  y = Genetic_Distance, color = Set),
              method = "lm", se = FALSE) +
  annotate("text", label = "All Samples", x = 0, y = 0.65, hjust = 0) +
  ylab(label = "Edward's Chord Genetic Distance") +
  xlab(label = "Distance Through Water") +
  xlim(c(0, max(IBD_All_Sets$Geographic_Distance)/1000)) +
  ylim(c(0, max(IBD_All_Sets$Genetic_Distance))) +
  theme(legend.position = "none") +
  scale_color_grey()

#### Unique MLGs ----

IBD_Plot_Unique_MLGs <- IBD_All_Sets %>% 
  filter(Distance_Type == "Distgen_noreps") %>% 
  ggplot() +
  geom_point(aes(x = Geographic_Distance/1000, 
                 y = Genetic_Distance, color = Set)) +
  geom_smooth(aes(x = Geographic_Distance/1000, 
                  y = Genetic_Distance, color = Set),
              method = "lm", se = FALSE) +
  annotate("text", label = "Unique MLGs", x = 0, y = 0.65, hjust = 0) +
  ylab(label = "Edward's Chord Genetic Distance") +
  xlab(label = "Distance Through Water") +
  xlim(c(0, max(IBD_All_Sets$Geographic_Distance)/1000)) +
  ylim(c(0, max(IBD_All_Sets$Genetic_Distance))) +
  theme(legend.position = "inside", 
        legend.position.inside = c(0.82, 0.3),
        legend.title = element_blank()) +
  scale_color_grey()

#### Combined Plot ----

IBD_Plot_Combined_by_Samples <-IBD_Plot_All_Samples/IBD_Plot_Unique_MLGs+
  patchwork::plot_layout(axis_titles = "collect")

ggsave(plot = IBD_Plot_Combined_by_Samples, here('figures','IBD_Plot_Combined_by_Samples.png'), width = 10, height = 15, units = "cm", dpi=400)

