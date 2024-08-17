if (!require("pacman")) install.packages("pacman")
pacman::p_load(diveRsity,
               here,
               tidyverse,
               reshape2)

#Load our custom functions - first check to see if
#it is loaded, and if not, load it.
if(!exists("summarize.divBasic.output", mode="function")) source(here::here("functions", "summarize.divBasic.output.r"))

if(!exists("partition.metric.pop.mean", mode="function")) source(here::here("functions", "summarize.divPart.metric.pop.mean.r"))

if(!exists("summarize.divPart.output", mode="function")) source(here::here("functions", "summarize.divPart.output.r"))

###Calculating basic popgen stats with the diveRsity package

####use genepop files created in the adegenet_create+dataframe file to calculate basic population genetic diversity stats.
#Depending on the size of your data set this step can take a while because it does the bootstraps for the allele rarefaction as the file is being read in.

#First create a list of all *.gen files you generated and want you run in diveRsity. These should be the noreps files only.

genepop.fileslist<-list.files("./genepop.files/", 
                               pattern="noreps", all.files=TRUE,
                               full.names=TRUE)

#get rid of the file path and the .gen suffix in the name list so you can use the data set name to name your new div. dataframes and later summary dataframes

genepop.filenames<-genepop.fileslist %>% 
  gsub('./genepop.files/', '',.) %>%
  gsub('.gen$', '',.)

#Read in list of geninds for this project and filter to keep only the ones that were used to make the genepop files you are using.  These will be used to join information about the sites back with the summarized pop gen data at the end of this script.

#Read in geninds for this project that were created and saved to an rds file using the adegenet_create_dataframe script.  

#To accommodate different names for the genind.rds file, you first use pattern matching to get the rds name, then read that in. You can only have one genind rds file in the output directory - the script will crash if you have more than one.

#After reading in the rds file, the geninds are filtered to keep only the ones that were used to make the genepop files you are using. These will be used to join information about the sites back with the summarized pop gen data at the end of this script. 

genind.list<-list.files("./output/", 
                        pattern=".genind.objects.rds", 
                        all.files=TRUE,
                        full.names=TRUE) %>% 
  readRDS(.) %>% 
  .[genepop.filenames] 


#Above, the . indicates the nascent list of genind objects that is piped after it is read in and you are using single brackets to select only elements where names are in the genepop.filenames object.

## Calculate within-population diversity statistics with divBasic from the diveRsity package.
divBasic.output<-genepop.fileslist %>%
  purrr::map(divBasic, outfile = NULL, #run divBasic on every data set
             gp = 3, bootstraps = 10, 
             HWEexact = TRUE) %>% 
  setNames(genepop.filenames) %>% #name the list item with dataset name
  purrr::map(summarize.divBasic.output) #use our function to summarize

################################################################
#############Calculate the partitions among populations#########
################################################################

#map through and run fastDivPar on all the .gen in your list, naming the output partition with the input file name.

#The use custom function to summarize the values output by fastdivPart(theta, G'st Hedrick, Gst, and Jost's D to get a mean pairwise difference of each site from all other sites.

DivPart.output<-genepop.fileslist %>% 
purrr::map(fastDivPart, outfile = NULL, pairwise = TRUE,
                      fst = TRUE, bs_locus = FALSE, bs_pairwise = TRUE,
                      boots = 5, plot = FALSE, para = FALSE) %>% 
setNames(genepop.filenames)%>% 
purrr::map(summarize.divPart.output)

#####################################################################
#############Join DivBasic and fastDivPart data################
#####################################################################

popgentables<-lapply(1:length(DivPart.output),function(x) merge(divBasic.output[[x]],DivPart.output[[x]],by = "Pop",  SIMPLIFY = FALSE) ) %>% 
setNames(names(DivPart.output))

#Add strata in from the respective genind files.
popgentables.with.strata<- lapply(1:length(popgentables),function(x) merge(distinct(genind.list[[x]]@strata),divBasic.output[[x]],by.x = "NewPop", by.y="Pop", SIMPLIFY = FALSE) )%>% 
  setNames(names(popgentables))

popgentables.with.strata %>% 
purrr::iwalk(~ write_csv(.x, here::here("output",paste0("diveRsity.popgentable.", .y, ".csv"))))
