##This script creates all the genind and genpop objects and genepop files you need for  population genetic analysis starting with the gendata list created in adegenet_create_dataframes.r script. 
create.popgen.objects.and.files <- function(gendata.list) {
  
  if (!dir.exists("./genepop.files")) {
        # create the "my_new_folder
    dir.create("genepop.files")
  }
  if (!dir.exists("./output")) {
    # create the "my_new_folder
    dir.create("output")
  }

    if(!exists("create.geninds", mode="function")) source(here::here("functions", "create.geninds.r"))
##############################################
##         CREATE GENIND OBJECTS            ##
##############################################

#The code below creates genind objects for .gendata dataframes in the specified list.  

#convert all our dataframes to genind objects
#using our custom function (create.geninds) that
#uses aedgenet's df2genind. You are using purrr:map to create
#genind.objects for each of the raw data sets you
#need to summarize.

#Map our custom create.geninds function over our list of gendata
#objects (make sure you sourced the function when
#you loaded packages).

list.of.genind.objects <- map(list.of.gendata.dataframes, ~create.geninds(.x))

#change the names to remove gendata and get just the base data set name.
names(list.of.genind.objects)<- sub(".gendata", "", names(list.of.genind.objects))

#save out the list of genind objects to a file you can use in the future.
saveRDS(list.of.genind.objects, file = here::here("output","list.of.genind.objects.rds"))

#If you want to create a genind object out of a single gendata dataframe, you can do so by using our function with any one of the dataframes in the gendata dataframes list as follows:
#single.example<-create.geninds(list.of.gendata.dataframes$noreps.lowerhudson.gendata)

##############################################
##         CREATE GENPOP OBJECTS            ##
##############################################

list.of.genpop.objects <- map(list.of.genind.objects, ~adegenet::genind2genpop(.x),  process.other=TRUE)

saveRDS(list.of.genpop.objects, file = here::here("output","list.of.genpop.objects.rds"))

#############################################################
######### Write GenePop files for each genind object ####
####################################################################

#You need to create genepop format files from your genind objects for use in diverSity. This is different from the genpop objects created a above.

#This code relies on functions that you got at install_github("romunov/zvau") above when you installed and loaded packages.

#Your genind object has @strata defined - this is done automagically by our function. 

#Note that you need to have at least two individuals in all populations for diveRsity to work.  Our dataframes are fine but you could run into problems with other data sets.

#This loop goes through all the genind objects that are stored in the list list.of.genind.objects and creates a genepop output file for each using the Zvau function.


for(i in 1:length(list.of.genind.objects)){
  writeGenPop(list.of.genind.objects[[i]], file=paste("./genepop.files/",(names(list.of.genind.objects[i])), ".gen", sep=""), comment="nocomment")
  
}
listList <- list(list.of.genind.objects=list.of.genind.objects, 
                 list.of.genpop.objects=list.of.genpop.objects)
return(listList)

}
#if you just need to write a file or two, here is example code
#writeGenPop(<<genind.name>>, file=here::here("genepop.files", "genind.name.gen", comment="")
#
#writeGenPop(list.of.genind.objects$hudbase.2014.2015.mainpops.genind, file="X:/My Drive/2020_Hudson/hudbase.2014.2015.mainpops.gen", comment="")
