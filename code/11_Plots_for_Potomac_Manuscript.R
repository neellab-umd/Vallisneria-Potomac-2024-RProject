if (!require("pacman")) install.packages("pacman")
pacman::p_load(here,
               magrittr,
               tidyverse,
               ggplot2,
               ggrepel,
               patchwork)

## Usage --

## This script is set up to use within the RProject Potomac.Manuscript.2024.R.Project.Rproj. It plots the results of several other scripts in ggplot. Those scripts need to have been run already so that their results are in the processed.data folder. All the necessary intermedicate data files exist in the .processed.data folder in the project but they can be created again with the provided scripts that live in the code folder.


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
## Population Level MLGs ####
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXX

## Plot mean and distribution of GD values in tidal versus non-tidal environments.
### Read in MLG data.

poppr.data <-readr::read_csv(here("processed.data", 
                           "poppr.site.level.output.allpotomac.csv")) %>% 
  dplyr::mutate(Tide=ifelse(Tide =="NonTidal","Nontidal", Tide))


#transform data set from wide to long to get the four variables into one graph

Three.GDs.DF <- poppr.data %>% 
  tidyr::pivot_longer(c("GD", "GD.SW", "GD.Simp"),
               names_to = "Genotypic.Diversity.Name", 
               values_to = "Genotypic.Diversity.Values") %>%
  dplyr::mutate(Genotypic.Diversity.Name = 
           case_match(Genotypic.Diversity.Name,
                                   "GD" ~ "All MLGs", 
                                   "GD.SW" ~"Effective Shannon's",
                                   "GD.Simp" ~ "Effective Simpson's")) %>%
  dplyr::mutate(Related.Name = factor(Genotypic.Diversity.Name, 
                               levels = c("All MLGs", 
                               "Effective Shannon's", 
                               "Effective Simpson's")))

ggplot(Three.GDs.DF, aes(x = factor(Tide), y = Genotypic.Diversity.Values, fill = Genotypic.Diversity.Name))+
  labs(x ="Tidal Regime", y = "GD") +
  scale_x_discrete(limits=c("Nontidal", "Tidal")) +
  theme_classic()+
  theme(axis.title.x = element_text(size=10),
        axis.text.x  = element_text(size=10),
        axis.title.y = element_text(size=10),
        axis.text.y  = element_text(size=10),
        legend.title = element_blank(),
        legend.text = element_text(size=9),
        legend.position = c(0.8, 0.17),
        legend.background = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(0.5, "cm"),
        legend.spacing.x = unit(0.1, 'cm'),
        legend.spacing.y = unit(0.2, 'cm'))+
  geom_violin(scale="width")+
  #ggforce::geom_sina(scale="width")+
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x), 
               position = position_dodge(.895),
               shape = 18, 
               size = .8,
               show.legend = FALSE,
               geom = "pointrange")+
  scale_fill_manual(values = c("#FFFFFF", "grey80", "grey50", "#FFFFFF", "grey80", "grey50"))

ggsave(here('figures','Three.GDs.bytide.png'), width = 15, height = 11.5, units = "cm", dpi=400)

#XXXXXXXXXXXXXXXXXXXXXXXXXXXX
##   Relatedness Graphs  ####
#XXXXXXXXXXXXXXXXXXXXXXXXXXXX

## Read in summary data and create factors to enforce correct plotting orders

allsamples <-readr::read_csv(here("processed.data","related","allsamples",
                                  "wang.rW.and.rA.allsamples.by.sample.csv")) 

noreps <-readr::read_csv(here("processed.data","related","noreps",
                              "wang.rW.and.rA.noreps.by.sample.csv")) 

related.data <- allsamples %>%
  dplyr::left_join(noreps, by = c("FromSamp","FromSite","Tide")) %>% 
  dplyr::mutate(Tide = factor(Tide))

## Transform data set from wide to long to get the four relatedness categories variables onto the x axis of one graph.

## Then change relatedness labels so they show up in the legend with more reasonable names than the original variable names and in the order I want.

wang.DF.all.and.noreps <- related.data %>% 
  tidyr::pivot_longer(c("Wang.allreps.rW",
                        "Wang.allreps.rA", 
                        "Wang.noreps.rW",
                        "Wang.noreps.rA"),
                      names_to = "Related.Name", 
                      values_to = "Related.Data") %>%
  dplyr::mutate(Related.Name = case_match(Related.Name,
                                          "Wang.allreps.rW" ~ "All Samples, Within Sites", 
                                          "Wang.allreps.rA" ~"All Samples, Among Sites", 
                                          "Wang.noreps.rW" ~"Unique MLGs, Within Sites",
                                          "Wang.noreps.rA" ~ "Unique MLGs, Among Sites")) %>%
  dplyr::mutate(Related.Name = factor(Related.Name, 
                                      levels = c("All Samples, Within Sites",
                                                 "All Samples, Among Sites", 
                                                 "Unique MLGs, Within Sites",
                                                 "Unique MLGs, Among Sites")))


ggplot(wang.DF.all.and.noreps, aes(x = factor(Tide), y = Related.Data, fill = Related.Name)) +
  labs(x ="Tidal Regime", y = "Wang's r") +
  scale_x_discrete(limits=c("NonTidal", "Tidal"),
                   labels = c("Nontidal", "Tidal")) +
  geom_violin(scale="width") + 
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x), 
               position = position_dodge(.895),
               shape = 18, 
               size = .8,
               show.legend = FALSE,
               geom = "pointrange",
               na.rm = TRUE) +
  theme_classic() +
  theme(axis.title.x = element_text(size=10),
        axis.text.x  = element_text(size=10),
        axis.title.y = element_text(size=10),
        axis.text.y  = element_text(size=10),
        legend.title = element_blank(),
        legend.text = element_text(size=9),
        legend.position = "inside",
        legend.position.inside = c(0.72, 0.84),
        legend.background = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(0.5, "cm"),
        legend.spacing.x = unit(0.1, 'cm'),
        legend.spacing.y = unit(0.2, 'cm')) +
  scale_fill_manual(values = c("#FFFFFF", "grey80", "grey50", "gray30","#FFFFFF", "grey80", "grey50", "gray30"))


ggsave(here('figures','Wang.within.and.amongpop.bytide.png'), width = 15, height = 11.5, units = "cm", dpi=400)


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
##      Correspondence Analysis      ####
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

ca.coords.noreps <- readr::read_csv("./processed.data/potomac.noreps.CA.csv")

## Plot of Correspondence Analysis results for axis 1 and 2 -no reps
ggplot(ca.coords.noreps, aes(x = CA.Dim.1, y = CA.Dim.2,label=OrderPop)) +
  geom_point(aes(color=Tide, shape=as.character(Year.Collected)), size=2.5,show.legend=FALSE)+
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  labs(x ="CA Axis 1, 30.2%", y = "CA Axis 2, 14.7%") +
  theme_classic()+
  theme(axis.title.x = element_text(size=9),
        axis.text.x  = element_text(size=9),
        axis.title.y = element_text(size=9))+
  scale_color_grey(start =0, end = .6)+
  geom_text_repel(show.legend=FALSE, size=3,min.segment.length = .41, max.overlaps = 14)

ggsave(here("figures",'CA.png'), width = 12, height = 12, units = "cm", dpi=400)

## Axes 3 & 4, not used in manuscript.

ggplot(ca.coords.noreps, aes(x = CA.Dim.3, y = CA.Dim.4,label=OrderPop)) +
  geom_point(aes(color=Tide, shape=as.character(Year.Collected)), size=2.5,show.legend=FALSE)+
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  labs(x ="CA Axis 3, 9.0%", y = "CA Axis 3, 6.7%") +
  theme_classic()+
  theme(axis.title.x = element_text(size=9),
        axis.text.x  = element_text(size=9),
        axis.title.y = element_text(size=9))+
  scale_color_grey(start =0, end = .6)+
  geom_text_repel(show.legend=FALSE, size=3,min.segment.length = .41, max.overlaps = 14)


