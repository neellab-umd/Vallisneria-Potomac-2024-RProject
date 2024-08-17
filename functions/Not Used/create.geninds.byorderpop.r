## adegenet Code to create genind and genpop objects##

create.geninds.by.order.pop <- function(gendata){
  # get the numbers for each column

  # create object with only your allele data to used by dfgenind
  gendata.alleles <- dplyr::select(gendata,aagx030:m16)
  
  strata_df <- gendata %>% 
    dplyr::select(OrderPop, NewPop, Tide) %>% 
    dplyr::mutate(OrderPop = as.factor(OrderPop),
           NewPop = as.factor(NewPop),
           Tide =as.factor(Tide)
           )
  
  # Read data matrix into a genind object including loci and population names
  # based on the column numbers for the loci as pulled by sclero
  gendind.ind <- df2genind(gendata.alleles,
                         ploidy=2, sep=":", 
                         type = "codom", 
                         ind.names=gendata$IDName,
                         pop=as.factor(gendata$OrderPop),
                         strata=strata_df
                         )
  
gendind.ind@other$latlong = dplyr::select(gendata,Longitude:Latitude)
gendind.ind@other$utms = dplyr::select(gendata,X:Y)

gendind.ind@other$CLONE.ID.2018<- gendata$Clone.ID.2018
gendind.ind@other$Tide<- as.factor(gendata$Tide)
gendind.ind@other$OrderPop<- as.factor(gendata$OrderPop)
gendind.ind@other$NewPop<- as.factor(gendata$NewPop)
gendind.ind@other$Tide_Pop<- paste0(gendata$Tide, sep="_",gendata$NewPop)
#strata(gendind.ind) <- data.frame(other(gendind.ind))

#

return(gendind.ind)
}