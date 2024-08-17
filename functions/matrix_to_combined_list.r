## When you convert the matrix of distances to a pairwise list, the some observations are only in the "to" column and other sites are only in the "from" column.  Those pairs for each of those elements will be left out of any summaries done by one direction. To get every site in our summary calculations, we have to put every pair in the list twice - once in the from:to direction and once in the to:from direction. Because you summarize by observations in only on one direction, you pick up every sample with no duplication. To do that we reverse the order of the samples and then rbind them together. 

## Given a starting matrix and the type of distance represented in that matrix, this function creates the combined pairwise list with all elements both in the from and to columns. The argument disttype is used to create the name for the column containing each distance

matrix_to_combined_list <- function(dist.matrix,disttype="Dist") {

  #Turn the matrix into columns and filter out self pairs
  sitedistlist<-setNames(reshape2::melt(as.matrix(dist.matrix)), c('To', 'From', 'Dist'))%>% 
    dplyr::rename(!!quo_name(disttype) := Dist) %>% 
    dplyr::filter(To!=From) #remove self pairs

  #pull the samples going from the To:From direction
  sitedistlist_to<-sitedistlist[c("To","From",disttype)]
  colnames(sitedistlist_to)<-c("To","From",disttype)
  
  #Reverse the order - put the From in the first column and the To in the second column
  sitedistlist_from<-sitedistlist[c("From","To",disttype)]
  colnames(sitedistlist_from)<-c("To","From",disttype)
  
  #Combine the from and to data sets using rbind
  sitedistlist_combined<-rbind(sitedistlist_to,sitedistlist_from)
  
  return(sitedistlist_combined)
}