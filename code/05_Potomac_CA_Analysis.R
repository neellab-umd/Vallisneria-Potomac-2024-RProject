## Author: Maile Neel
## 2024

## Usage ---

## The script is designed to work within the RProject Potomac.Manuscript.2024.R.Project.Rproj

#It implements Correspondence Analysis using a genpop objects created in the 01.PR_Manuscript_adegenet_Create_Dataframes.r script.  If that script has been run, the genpop objects will be in a list called list.of.genpop.objects which is saved to an rds file called list.of.genpop.objects.rds in the ./processed.data folder. Make sure all those objects first before proceeding here.

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Install and Load Packages ####
##XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

if (!require("pacman")) install.packages("pacman")
pacman::p_load(here,
               tidyverse,
               adegenet,
               factoextra,
               FactoMineR)

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#      Read in Data      ####
##XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

## If the list.of.genpop.objects.rds file is not in the processed.data folder, run the 01_Create_adegenet_geninds.R script first to generate the object.

list.of.genpop.objects <- readRDS(here("processed.data","list.of.genpop.objects.rds"))

##Extract additional site variables from the genpop that will be used at the end to join back to the CA results. Getting all our site variables back out of the genpop object is a pain because genind2genpop() does not collapse categorical variables in the "other" slot down to one per population.  We have to do it manually.

site.labels <-  cbind(as.character(list.of.genpop.objects$genpop.allpotomac@other$OrderPop),
                      as.character(list.of.genpop.objects$genpop.allpotomac@other$NewPop),
                      as.character(list.of.genpop.objects$genpop.allpotomac@other$Tide),
                      as.character(list.of.genpop.objects$genpop.allpotomac@other$Year.Collected)) %>% 
  as.data.frame() %>% 
  dplyr::rename(OrderPop = V1, NewPop = V2, Tide = V3, Year.Collected = V4) %>% 
  dplyr::group_by(NewPop) %>% 
  dplyr::summarize(OrderPop = first(OrderPop),
                   NewPop = first(NewPop),
                   Tide = first(Tide),
                   Year.Collected = first(Year.Collected)) 

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Correspondence Analysis Only Unique MLGs   ####
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

# Use FactoMineR to do Correspondence analysis on no reps data set.

## http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/113-ca-correspondence-analysis-in-r-essentials/#r-code-to-compute-ca

res.ca.noreps <- FactoMineR::CA(list.of.genpop.objects$genpop.noreps.potomac@tab, graph = FALSE)

summary(res.ca.noreps)

ca.coords.noreps <- factoextra::get_ca_row(res.ca.noreps)$coord %>% 
  as.data.frame() %>% 
  rownames_to_column("NewPop") %>% 
  rename("CA.Dim.1"= "Dim 1", "CA.Dim.2"= "Dim 2", 
         "CA.Dim.3"= "Dim 3", "CA.Dim.4"= "Dim 4",
         "CA.Dim.5"= "Dim 5") %>% 
  left_join(site.labels, by = "NewPop")

#++++++++++++++++++++++++
##   Write out CA scores    ####
#++++++++++++++++++++++++

## Write out CA scores for later plotting with ggplot2.

write_csv(ca.coords.noreps,"./processed.data/potomac.noreps.CA.csv")

## Code for plotting the Correspondence Analysis results is in the script 11_Plots_for_Potomac_Manuscript.R

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Correspondence Analysis All Samples  ####
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

## This was run for exploratory purposes but then was not used in the manuscript. Only the CA based on unique MLGs was used.

res.ca.allsamples <- FactoMineR::CA(list.of.genpop.objects$genpop.allpotomac@tab, graph = FALSE)

ca.coords.allsamples <- factoextra::get_ca_row(res.ca.allsamples)$coord %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("NewPop") %>% 
  dplyr::rename("CA.Dim.1"= "Dim 1", "CA.Dim.2"= "Dim 2", 
                "CA.Dim.3"= "Dim 3", "CA.Dim.4"= "Dim 4",
                "CA.Dim.5"= "Dim 5") %>% 
  dplyr::left_join(site.labels, by = "NewPop")

## Plot of Correspondence Analysis results for axis 1 and 2 - all samples

ggplot(ca.coords.allsamples, aes(x = CA.Dim.1, y = CA.Dim.2,label=OrderPop)) +
  geom_point(aes(color=Tide, shape=as.character(Year.Collected)), size=2.5,show.legend=FALSE)+
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  labs(x ="CA Axis 1, 32.7%", y = "CA Axis 2, 21.4%") +
  theme_classic()+
  theme(axis.title.x = element_text(size=9),
        axis.text.x  = element_text(size=9),
        axis.title.y = element_text(size=9))+
  scale_color_grey(start =0, end = .6)+
  ggrepel::geom_text_repel(show.legend=FALSE, 
                           size=3,min.segment.length = .41, max.overlaps = 14)

ggsave(here("figures",'CA.png'), width = 12, height = 12, units = "cm", dpi=400)

ggplot(ca.coords.allsamples, aes(x = CA.Dim.3, y = CA.Dim.4,label=OrderPop)) +
  geom_point(aes(color=Tide, shape=as.character(Year.Collected)), size=2.5,show.legend=FALSE)+
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  labs(x ="CA Axis 3, 11.3%", y = "CA Axis 3, 7.4%") +
  theme_classic()+
  theme(axis.title.x = element_text(size=9),
        axis.text.x  = element_text(size=9),
        axis.title.y = element_text(size=9))+
  scale_color_grey(start =0, end = .6)+
  geom_text_repel(show.legend=FALSE, size=3,min.segment.length = .41, max.overlaps = 14)


