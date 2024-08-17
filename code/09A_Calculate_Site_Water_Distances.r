## Maile Neel
## 2024

## Usage ----


## This script is set up to use within the  RProject Potomac.Manuscript.2024.R.Project.Rproj. It calculates cost distances through water for the sites in the Potomac River using a transition matrix created in the script 07_Create_gdistance_cost_surface.r The transition matrix is in the rds file tr1c.Potomac.Transition.10m.rds in the data/cost_surfaces directory.  If the transition matrix does not exist in that directory, that script needs to be run to create the transition matrix prior to running this script. 
 
## The script also calculates Euclidean distances among sites. To reduce the effects of error imposed by raster-based cost distances underestimating actual distance when sites are close to one another, we take the larger of the cost or Euclidean distance between sites as the final distance.

## Results are output to four files. The following three RDS files of distance matrices (one for the whole river, and one each for tidal and non-tidal regions) are used in the script 10_Potomac_IBD_Among_Sites.r.
## Full_PR_Site_Cost_and_Euclidean_Distances_combined_Matrix
## Nontidal_PR_Site_Cost_and_Euclidean_Distances_combined_Matrix
## Tidal_PR_Site_Cost_and_Euclidean_Distances_combined_Matrix

##  The csv called PR_Site_Cost_and_Euclidean_Distances_combined.csv is used in the script 09B_Summarize_Site_Water_Distances.r. All are in the processed.data folder.  If those files are present, this script does not been to be run again.

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Install and Load Packages and Source Functions ----
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,
               here, 
               gdistance,
               terra,
               sf, 
               ggthemes,
               reshape2,
               update = FALSE)

if(!exists("matrix_to_combined_list", 
           mode="function")) source(here::here("functions",
                                               "matrix_to_combined_list.r"))

## Steps 
## 1. Read in point data.
## 2. Read in Potomac River shapefile for plotting purposes.
## 3. Read in Potomac River transition matrix.
## 4. Calculate cost distances.
## 5. Process distance matrix.


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Read In Site Centroid Point Data ----
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

## In this case the points are the sample locations used in the Potomac River manuscript. 
## To work with the genetic matrices when calculating IBD, the sites need to be alphabetical, so we are arranging by NewPop.

Potomac.sites <- read_csv(here::here("data","PR.Pop.Centroids.UTMs.csv")) %>% 
 dplyr::arrange(NewPop) %>% 
  dplyr::select(-c(Longitude,Latitude)) #Get rid of these coordinates
  
  
#We  project the points into UTM with the NAD83 datum (epsg:3725) so that the points are in the same projection as our watermask. We use the Cartesian UTM system to get accurate distance measurements.

Potomac.sites.sf.utm <- Potomac.sites %>% 
  column_to_rownames("NewPop") %>% 
  sf::st_as_sf(                   #convert the points to an sf object
  coords=c("X", 
           "Y"),
  agr = "constant",
  crs = 26918,               #project to UTM NAD83
  stringsAsFactors = FALSE,
  remove = TRUE
) 

Potomac.tidal.sites.sf.utm <- Potomac.sites.sf.utm %>% 
  dplyr::filter(Tide == "Tidal")

Potomac.nontidal.sites.sf.utm <- Potomac.sites.sf.utm %>% 
  dplyr::filter(Tide == "Nontidal")

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Read in shapefile for visualization----
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

## You only need to do this if you want to see the point locations relative to the River.  It is not used in the calculations.

## Read in Chesapeake watermask shapefile and crop it to a bounding box that includes all Potomac sample sites plus a bit of space to include all the water between them.

Potomac_watermask.polygon.utm <-  
  sf::st_read(here::here("shapefiles",
                         "CB.PR.Digitized.MN.2024.shp")) %>% 
  sf::st_transform(
    crs = 26918 #transform to utm NAD83 
  ) %>% 
  dplyr::filter(Water == "Potomac")

## Plot shapefile and points for a quick check ----

ggplot() + 
  geom_sf(data = Potomac_watermask.polygon.utm, 
          size = 1, 
          color = "blue", 
          fill = "cyan1") + 
  ggtitle("Watermask with Potomac Sites") + 
  geom_sf(data = Potomac.sites.sf.utm,
          size = 2, 
          color = "magenta")+
  coord_sf()+
  theme_map()

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#  Load transition matrix -----
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#If it is not loaded already, read in the resistance surface you created in the PR_gdistance_surface_code. You will want to use the corrected version (indicated by the c after tr1 in the name). The projection correction is needed for accuracy.

tr1c_Potomac_10 <- readRDS(here::here("cost_surfaces","tr1c.Potomac.Transition.10m.rds"))

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Calculate Euclidean Distance Among All Site Pairs ----
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

## Create a least cost distance matrix between points use the coordinates of point locations read in above

## Cost is defined as the reciprocal of the values in the transition matrix

## the as(., "Spatial") function converts an sf object to the SpatialPointsDataFrame that is needed by gdistance's costDistance function

Site.rescost.dist10 <- gdistance::costDistance(tr1c_Potomac_10, as(Potomac.sites.sf.utm, "Spatial"))

Site.rescost.dist10.tidal <- gdistance::costDistance(tr1c_Potomac_10, as(Potomac.tidal.sites.sf.utm, "Spatial"))

Site.rescost.dist10.nontidal <- gdistance::costDistance(tr1c_Potomac_10, as(Potomac.nontidal.sites.sf.utm, "Spatial"))

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Calculate Euclidean Distance Among All Site Pairs ----
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

Site.Euclid.dist <- Potomac.sites %>% 
    column_to_rownames("NewPop") %>% #by moving the population name label to a rowname, dist will use it for column and row names
  dplyr::select(X,Y) %>% 
  dist(method = "euclidean", diag =TRUE, upper = TRUE)

  Site.Euclid.dist.tidal <- Potomac.sites %>% 
    dplyr::filter(Tide == "Tidal") %>% 
    column_to_rownames("NewPop") %>% #by moving the population name label to a rowname, dist will use it for column and row names
    dplyr::select(X,Y) %>% 
    dist(method = "euclidean", diag =TRUE, upper = TRUE)
  
    Site.Euclid.dist.nontidal <- Potomac.sites %>% 
      dplyr::filter(Tide == "Nontidal") %>% 
    column_to_rownames("NewPop") %>% #by moving the population name label to a rowname, dist will use it for column and row names
      dplyr::select(X,Y) %>% 
      dist(method = "euclidean", diag =TRUE, upper = TRUE)
  

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#  Get the larger value between cost and Euclidean -----
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

Site.Final.dist <- pmax(Site.Euclid.dist,Site.rescost.dist10)

Site.Final.dist.tidal <- pmax(Site.Euclid.dist.tidal,Site.rescost.dist10.tidal )
  
Site.Final.dist.nontidal <- pmax(Site.Euclid.dist.nontidal,Site.rescost.dist10.nontidal )

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#  Write Distance Matrices To RDS Files-----
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

saveRDS(Site.Final.dist, here("processed.data", "Full_PR_Site_Cost_and_Euclidean_Distances_combined_Matrix_Ordered.by.rds"))

saveRDS(Site.Final.dist.tidal, here("processed.data", "Tidal_PR_Site_Cost_and_Euclidean_Distances_combined_Matrix_Ordered.by.rds"))

saveRDS(Site.Final.dist.nontidal, here("processed.data", "Nontidal_PR_Site_Cost_and_Euclidean_Distances_combined_Matrix_Ordered.by.rds"))


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
### Turn Final Distance Matrix into pairwise list of distances ----
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#Use custom function matrix_to_combined_list() to 1) turn the Site.Final.dist matrix into columns and 2) Give the type of distance (i.e., Geog for geographic, or Water for over water only, or Final) 3) filter out self pairs and make it so every pair of elements in the distance matrix shows up in the list twice - once in the from:to direction and once in the to:from direction. Because you summarize only on one direction, you pick up every sample with no duplication. 

## Also add a column that flags samples as being from the same versus different tidal environment. The flags are used in summarizing in the script 09B_Summarize_Site_Water_Distances.r. They are probably not necessary but it made later summaries much easier for me. 


site.fnal.distlist_combined <- 
  matrix_to_combined_list(Site.Final.dist, disttype = "FinalDist") %>%
  dplyr::rename(ToPop = To, FromPop = From) %>% 
  left_join(Potomac.sites, by = c("ToPop" = "NewPop")) %>% 
  rename_with(.cols = OrderPop:Y, ~paste0("To_",. )) %>% 
  left_join(Potomac.sites, by = c("FromPop" = "NewPop")) %>%
  rename_with(.cols = OrderPop:Y, ~paste0("From_",. ))%>% 
  mutate(FinalDist = round(FinalDist/1000,1)) %>%  #convert to km
  mutate(SameTide = ifelse(From_Tide == To_Tide, "Same","Diff")) 


dim(site.rescost.distlist_combined)
head(site.rescost.distlist_combined)

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#  Write Distances To csv File-----
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

write_csv(site.rescost.distlist_combined, here("outputs", "PR_Site_Cost_and_Euclidean_Distances_combined_List.csv"))
