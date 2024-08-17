#This function is used to iterate over pairwise metrics in the list that is output from fastDivPart in diveRsity and converts the matrices to pairwise lists, and then calculates the mean value of each pairwise metric.  
#
#It returns a list that contains one dataframe for each metric calculated.
#
#It is used within another custom function summarize.divPart.output() that further cleans up the values to create a single summary table.

#Usage is
#partition.metric.pop.mean(input.object)


partition.metric.pop.mean<-function(metric) {
  
  distlist<-as_tibble(metric, 
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
    dplyr::summarise(metric = mean(value))
  
  #return(metricbypop)

}