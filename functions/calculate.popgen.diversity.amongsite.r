#This function is used to iterate over pairwise metrics in the list that is output from mmod and hierfstat and converts the matrices to pairwise lists of sites, and then calculates the mean value of each pairwise metric for each site.  
#
#It returns a dataframe of the mean values that contains one dataframe for each metric calculated.
#
#It is used within another custom function that further cleans up and aggregates the among-site diversity statistics to create a single summary table.

#Usage is
#calculate.popgen.diversity.amongsite(input.object)


calculate.popgen.diversity.amongsite<-function(outlist, metric,...) {
  
  distlist<-as_tibble(as.matrix(outlist), 
                      rownames = 'Pop') %>% 
    tidyr::pivot_longer(names_to = "PopFrom", 
                        values_to = "value",
                        -Pop, 
                        values_drop_na = TRUE) %>% 
    dplyr::mutate_all(~sub(",", "", .)) %>% 
    dplyr::mutate(value = as.numeric(value))
  
  distlist_Reverse<-distlist[,c(2,1,3)]
  
  colnames(distlist_Reverse)<-c("Pop", "PopFrom", "value")
  
  fullmetric<-rbind(distlist,distlist_Reverse)
  
  metricbypop<-as.data.frame(fullmetric) %>%
    dplyr::group_by(Pop) %>%
    dplyr::summarise(!!metric := round(mean(value),3))
  
  #return(metricbypop)
  
}
