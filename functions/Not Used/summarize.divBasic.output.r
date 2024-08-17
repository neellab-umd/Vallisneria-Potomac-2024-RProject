#This function extracts the overall basic popgen data from that comes from divBasic and puts it into a table.
# 

#Get just the overall statistics into a list in which each population is an element of the list

  summarize.divBasic.output <- function(diveRsity.data) {
    overall.list<-(lapply(diveRsity.data$mainTab,"[", 1:10,c("overall"),drop=FALSE))
    
    #flatten/unlist the elements of overall.list within each population
    overall.flattened<-lapply(overall.list, unlist)
    
    #convert to a dataframe
    df.overall<-as.data.frame(overall.flattened)
    
    #transpose the dataframe so the statistics are in columns
    popgentable<-as.data.frame(t(df.overall))
    
    colnames(popgentable)<-c("N", "A","%","Ar","Ho","He","HWE","Fis","Fis_Low","Fis_High")
    
    #convert row names to Pop variable and put it in the first column. #get rid of the underscore and number genepop format uses in the population name, and create a more sensible Fis CI variable.
    popgentable<-tibble::rownames_to_column(popgentable, 
                                            var = "Pop") %>% 
      separate(Pop, c("Pop"), sep = "_", extra='drop') %>% 
      mutate(Fis_CI = paste(Fis_Low,Fis_High, sep = " - ")) %>% 
      dplyr::select(-Fis_Low,-Fis_High)
  
    
  #****assign popgentable to a dataframe with a more meaningful name reflecting the input file
    #using the name you gave up above when you selected the dataframe
    assign(paste0("popgentable.", "Croton", sep=""),popgentable )
  }