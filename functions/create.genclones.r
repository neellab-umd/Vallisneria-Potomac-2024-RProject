##adegenet Code to create genind and genpop objects##

#Creates a genclone object from previously created genind objects using the CLONE.ID.2018 slot as the mlg identifiers. Defining the mlgs ahead of time (mlgclass = TRUE) prevents poppr from calling new mlgs each time it is run. This allows us to keep our naming scheme.  

create.genclones <- function(gendind.ind){
  #get the numbers for each column

genclone.gc<-as.genclone(gendind.ind, mlg=gendind.ind@other$CLONE.ID.2018, mlgclass = TRUE)


return(genclone.gc)
}