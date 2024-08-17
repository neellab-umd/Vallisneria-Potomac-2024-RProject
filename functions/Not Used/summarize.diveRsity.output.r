#This function extracts the basic popgen data from the list output by diveRsity after running divBasic. Basic stats has a slightly different format - the list label for that is main_tab (rather than mainTab) and there are only 9 variables rather than the 10 here.  


#Get just the overall statistics into a list in which 
#each population is an element of the list
summarize.diveRsity.output<- function(diveRsity.data) {
  overall.list<-(lapply(diveRsity.data$mainTab,"[",1:10,c("overall"),drop=FALSE))
  
  #flatten/unlist the elements of overall.list within each population
  overall.flattened<-lapply(overall.list, unlist)
  
  #convert to a dataframe
  df.overall<-as.data.frame(overall.flattened)
  
  #transpose the dataframe so the statistics are in columns
  popgentable<-as.data.frame(t(df.overall))
  popgentable<-tibble::rownames_to_column(popgentable)
  colnames(popgentable)<-c("Pop","N", "A","Perc_A","Ar", "Ho","He","HWE","Fis","Fis_Low","Fis_High")
  
    #get rid of the underscore and number genepop format uses in the population name you can't just use an underscore because it eliminates the _A and _B in BB and PTR names. The separators are specific to the Potomac data set given how genepop and DiverSity name populations.
  
  popgentable2<-popgentable %>% 
    separate(Pop, "Pop", sep = "_0", extra='drop') %>% 
    separate(Pop, "Pop", sep = "_1", extra='drop') %>% 
    separate(Pop, "Pop", sep = "_P", extra='drop')
  
  
  return(popgentable2)
}
