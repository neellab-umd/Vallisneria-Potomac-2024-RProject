## Author: Maile Neel

## Usage ----

## This script generates both Euclidean and cost distances among samples from the Potomac River. Both types of distances are used because resistance distances are not good for individuals within sites or more generally for points that are close to one another.  They are most problematic when the inter-sample distance is less than the minimum cell size in the resistance surface (10 m in our case). 

## Samples that are within <10m will be assigned a cost distance of 0 by the gdistance::costDistance() function. Samples that are between 10m and 20m apart will be assigned a distance of 10m. This assignment is the result of the fact that gdistance uses a centroid to centroid measurement of cell to cell movement. Anything less than 10 m is considered to be within one cell, all measurements are in effect 10 m less than they "should" be. The effect is most pronounced at short distances, but all cost distances will still show an effect of resistance surface cell size if there is no barrier to direct movement. At longer distances the cost distance becomes longer and more accurate due to bends in the river.

## For the most part, samples that are close to one another are in the same site. To limit the effect of error in cost distances, we use Euclidean distance for samples from the same site.

## There are also samples that are at different sites because they were sampled in different years, but samples are actually very close together.  In those cases we take the larger of the two distances because where Euclidean distance is greater than cost distance, it is more accurate and thus, it is retained in the final combined euclidean+cost distance dataframe. 

## The output of this script is the csv of all pairwise distances among samples with all the categories needed for summary, analysis and plotting that is done in 08B_Summarize_Sample_Distances.R.

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Install and Load Packages and Source Functions ----
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,
               here,
               reshape2,  ## used in matrix_to_combined_list() function 
               gdistance,
               terra,
               sf,
               update = FALSE
               )


if(!exists("matrix_to_combined_list", 
           mode="function")) source(here::here("functions",
                                               "matrix_to_combined_list.r"))


## To make this easily generalizeable to different subsets of our Vallisneria genetic data, the input dataframe is assigned to a dataframe called "maindf".  That dataframe has to have UTM coordinates in meters in columns labeled X and Y, the MLG ID in the variable called Clone.ID.2018, and the sample name labeled in IDName.  The code and custom function can be altered to allow for different variable names.

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Read In Sample Point Data ----
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

## Read in csv file for the Potomac samples then create a variable we will use for pulling different sets of comparisons.  

allpotomac <- read.csv(here::here("data",
                                       "allpotomac.microsatellite.data.csv")) %>% 
  dplyr::filter(Clone.ID.2018 != "NA") %>% 
  dplyr::select(IDName,NewPop,OrderPop,Clone.ID.2018,Longitude,Latitude, X, Y, Tide) %>% 
  mutate(poptideindmlg = paste0(as.character(Tide), ".",
                               as.character(NewPop),".",
                               as.character(IDName),".",
                               as.character(Clone.ID.2018))) %>% 
  rename(MLG = Clone.ID.2018,
         Pop = NewPop)


## assign allpotomac to maindf.  This step allows you to change the input dataset without replacing the name throughout the code so you do not have to repeat everything for a different data set.

maindf_nomissing<-allpotomac

## name the dataset you are using so the final output dataframes and csvs are named
## appropriately.

datasetname<-"allpotomac"

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Calculate Euclidean Distance Among All Samples ----
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

## Create a distance matrix of pairwise Euclidean distances among all individuals in meters based on UTM coordinates. To add column and row headings you need to use as.matrix to turns the distance matrix into a real matrix.

sampledist<- as.matrix(dist(maindf_nomissing[c("X","Y")], 
                            method = "euclidean", 
                            diag = FALSE, 
                            upper=FALSE), labels=TRUE)

## add the sample IDName to identify the
## samples that are linked to each matrix cell

colnames(sampledist) <- rownames(sampledist) <- maindf_nomissing$IDName

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
### Turn matrix into pairwise list of distances ----
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

## Use the custom function matrix_to_combined_list() to 1) turn the matrix into columns and 2) Give the type of distance (i.e., Geog for geographic, or Water for over water only) 3) filter out self pairs and make it so every pair of elements in the distance matrix shows up in the list twice - once in the from:to direction and once in the to:from direction. Because you summarize only on one direction, you pick up every sample with no duplication. 

## We also add column that flag samples as being from the same versus different MLG, Site, and tidal regime. These flags are used in summarizing below. They are probably not necessary but it made filtering and grouping later much easier for me. 

sample.geog.distlist_combined <- 
  matrix_to_combined_list(sampledist,disttype = "GeogDist") %>%
  dplyr::rename(ToSample = To, FromSample = From) %>% 
  dplyr::left_join(dplyr::select(maindf_nomissing, 
                          IDName,Pop,OrderPop,MLG,Tide), 
            by = c("ToSample" = "IDName")) %>% 
  dplyr::rename_with(.cols = Pop:Tide, ~paste0("To",. )) %>% 
  dplyr::left_join(dplyr::select(maindf_nomissing, 
                          IDName,Pop,OrderPop,MLG,Tide), 
            by = c("FromSample" = "IDName")) %>% 
  dplyr::rename_with(.cols = Pop:Tide, ~paste0("From",. ))%>% 
  dplyr::mutate(GeogDist=round(GeogDist,2),
                SameTide = 
                  ifelse(FromTide==ToTide, 
                         "Same","Diff"),
                SameSite = 
                  ifelse(FromPop == ToPop,
                         "Same","Diff"),
                SameMLG = 
                  ifelse(FromMLG == 
                           ToMLG,"Same","Diff"))
  
## Check the object

sample.geog.distlist_combined
dim(sample.geog.distlist_combined)
head(sample.geog.distlist_combined)
names(sample.geog.distlist_combined)

## This dataframe is now ready for joining with the dataframe of cost distances ## that is generated below

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#  Calculate Cost Distance Among Samples ######## 
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

### Read in the transition matrix ----

## If it is not loaded already, read in the resistance surface you created in the 07_Create_gdistance_cost_surface.r script. When this matrix was created I used a custom function to ensure that all the points were within water - if they are not, the cost distance function will fail. 

## You will want to use the corrected version (indicated by the c after tr1 in the name). The projection correction is needed for accuracy. 

tr1c.potomac10 <- readRDS(here::here("cost_surfaces","tr1c.Potomac.Transition.10m.rds"))

### Convert the points dataframe into a spatial sf object ----

## Convert the points to an sf object and project them into UTM with the NAD83 datum (epsg:26918) so that the points are in the same projection as the transition matrix.  We use the Cartesian UTM system to get accurate distance measurements.

maindf_nomissing.individuals.sf.utm <- maindf_nomissing %>% 
  sf::st_as_sf( 
    coords=c("X", 
             "Y"),
    agr = "constant",
    crs = 26918, ## project to UTM NAD83
    stringsAsFactors = FALSE,
    remove = TRUE
  ) 

## check the structure of the object

maindf_nomissing.individuals.sf.utm

### Create a least cost distance matrix between points ----

## Cost is defined as the reciprocal of the values in the transition matrix.

## the as(., "Spatial") function converts the sf object to the SpatialPointsDataFrame that is needed by gdistance's costDistance function

allpotomac_individual_rescost_Dist10 <- 
  costDistance(tr1c.potomac10, 
               as(maindf_nomissing.individuals.sf.utm,
                  "Spatial")) %>% 
  as.matrix(., labels=TRUE)


## Convert the distance matrix to a regular matrix and add sample IDs to the matrix row and column names.

colnames(allpotomac_individual_rescost_Dist10) <- rownames(allpotomac_individual_rescost_Dist10) <-maindf_nomissing.individuals.sf.utm$IDName

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
### Turn matrix into pairwise list of distances ----
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

## Use custom function to turn matrix into a pairwise "list" of samples 

sample.rescost.distlist_combined_Dist10 <- 
  matrix_to_combined_list(allpotomac_individual_rescost_Dist10,
                          disttype = "WaterDist") %>%
  dplyr::rename(ToSample = To, FromSample = From) 
  
dim(sample.rescost.distlist_combined_Dist10)
head(sample.rescost.distlist_combined_Dist10)

## Double check for points outside of water - i.e., distance of Inf.  Ideally this yields n() 0 given we modified the water mask to encompass all sample points. If there are any points with Inf distance, you need to go back to making sure all points are within water and that all the cells they are in are connected to the other water cells.

filter(sample.rescost.distlist_combined_Dist10, WaterDist=="Inf") %>% summarize(n())

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#   Join Euclidean and Cost Distance  ####
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

## Join the dataframes and furuthr process the columns to prep for analysis.

## For comparisons that include multiple sites, you will want cost distances. However, due to error involved in those distances when samples are close together, you need to use Euclidean distance. The final mutate in this pipe creates a variable that takes the cost distance if the sites are different and the Euclidean distance if the sites are the same. Croton needs to be managed manually because the two years are coded as different sites, but you need to use Euclidean distance. 

Ind.Euclidean.and.Cost.Dist.List<- sample.geog.distlist_combined %>% 
  relocate(GeogDist, .after = last_col()) %>% 
  bind_cols(WaterDist= sample.rescost.distlist_combined_Dist10$WaterDist) %>% 
  mutate(FinalDist = case_when(
    SameSite >= "Same" ~ GeogDist,## use Euclidean Distance within in sites
    str_detect(FromPop, "BWK") & str_detect(ToPop, "BWK") ~ GeogDist,
    str_detect(FromPop, "HCK") & str_detect(ToPop, "HCK") ~ GeogDist,
    str_detect(FromPop, "WSP") & str_detect(ToPop, "WSP") ~ GeogDist,
    TRUE ~ WaterDist ## Use cost distance for all other values 
  ))                ## not covered by the conditions above.

head(Ind.Euclidean.and.Cost.Dist.List)

Ind.Euclidean.and.Cost.Dist.List %>% 
  filter(FromPop != ToPop & WaterDist>FinalDist) %>% 
  filter(grepl("BWK", FromPop)) %>% 
  dplyr::select(-c(ToOrderPop:ToTide, FromOrderPop:SameMLG)) %>% 
  arrange(FinalDist)


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#   Save Distances to csv  ####
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

## This is a very large csv that cannot be fully loaded in Excel.

write_csv(Ind.Euclidean.and.Cost.Dist.List, file=here("processed.data", paste0(datasetname, ".euclidean.and.cost.distances.among.samples.csv")))

## This csv file is summarized in 8B_Summarize_Sample_Distances.R