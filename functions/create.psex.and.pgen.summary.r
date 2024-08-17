# Maile Neel 10/21/2021

## This  script creates the function called pgen.and.sex that is used to calculate, shockingly, pgen and psex for a given genind object.

#To run it, first source the function and then call pgen.and.sex(gid) where gid is your genind name that you have already loaded. The code where it is run for the Potomac manuscript is in the script 03_psex_poppr.R in the Potomac.Manuscript.RProject.

pgen.and.sex <- function(gid) {
  
  ##The pgen calculation is the probability of seeing a given multilocus genotype based on the allele frequencies for one rep of each MLG in each stratum indicated in the pop slot.  
  
  #Homozygote probabilities are p^2 and heterozygotes are 2pq.
  
  #+++++++++++++++++++++++++++++++++++
  #######  Calculate pgen   ##########
  #+++++++++++++++++++++++++++++++++++
  
  # Calculate the overall probability of the MLG by taking the product of the log probabilities of each locus and then taking the exponent.
  
  ## first within each population
  pgen.pop.exp <- pgen(gid, 
                       by_pop = TRUE) %>% # calculate pgen matrix
    rowSums(na.rm = TRUE) %>% # take the sum of each row
    exp() %>% # take the exponent of that sym
    as.data.frame() %>% 
    dplyr::as_tibble(rownames = "sample")
  colnames(pgen.pop.exp) <- c("sample","Pgen.within.pop")
  
  ## then overall  
  pgen.all.exp <- pgen(gid, 
                       by_pop = FALSE) %>% # calculate pgen matrix
    rowSums(na.rm = TRUE) %>% # take the sum of each row
    exp() %>% # take the exponent of the results
    as.data.frame() %>% 
    dplyr::as_tibble(rownames = "sample")
  colnames(pgen.all.exp) <- c("sample","Pgen")
  
  pgen.both.exp <-pgen.all.exp %>% 
    left_join(pgen.pop.exp, by="sample")
  
  #+++++++++++++++++++++++++++++++++++++++++
  #######  Calculate psex single  ##########
  #+++++++++++++++++++++++++++++++++++++++++
  
  #Calculate the probability of a second sample of the same MLG being the result of sexual reproduction.
  ## First by population
  psex.second.bypop <- psex(gid, 
                            by_pop = TRUE,
                            method = "single") %>% 
    as.data.frame(.) %>% 
    dplyr::as_tibble(rownames = "sample")
  colnames(psex.second.bypop) <- c("sample","Psex.second.within.pop")
  
  ## Then overall
  
  psex.second <- psex(gid, 
                      by_pop = FALSE,
                      method = "single") %>% 
    as.data.frame(.) %>% 
    dplyr::as_tibble(rownames = "sample")
  colnames(psex.second) <- c("sample","Psex.second")
  
  psex.second.both <-psex.second %>% 
    left_join(psex.second.bypop, by="sample")
  
  #+++++++++++++++++++++++++++++++++++++++++++
  #######  Calculate psex multiple  ##########
  #+++++++++++++++++++++++++++++++++++++++++++
  
  #Calculate the probability of the MLG being seen as many times as it was in the number of samples taken if the observation was due to sexual reproduction.  
  
  #Calculate for each population separately - for this calculation the first encounter in each population is a probability at or near 1. So this part of the function only tells you about multiple encounters in that population.
  
  # Calculated by population
  psex.multiple.bypop <- psex(gid, 
                              by_pop = TRUE, 
                              method = "multiple"
  ) %>% 
    as.data.frame(.) %>% 
    dplyr::as_tibble(rownames = "sample")
  colnames(psex.multiple.bypop) <- c("sample","Psex.multiple.within.pop")
  
  # Then calculated regardless of population
  psex.multiple <- psex(gid, 
                        by_pop = FALSE, 
                        method = "multiple"
  )  %>% 
    as.data.frame(.) %>% 
    dplyr::as_tibble(rownames = "sample")
  colnames(psex.multiple) <- c("sample","Psex.multiple")
  
  
  psex.multiple.both <-psex.multiple %>% 
    left_join(psex.multiple.bypop, by="sample")
  
  
  #++++++++++++++++++++++++++++++++++++++++++++++++++
  ####### Join dataframes for final output ##########
  #++++++++++++++++++++++++++++++++++++++++++++++++++
  
  #Create a final data frame that includes pgen, psex.second. Psex.multiple, the mlg ID, number of times the mlg was seen, total number of samples and mlgs.
  
  #First get the numbers of mlgs and samples in each population for merging into the final dataframe.
  NumMLGs <- rowSums(mlg.table(gid, plot = FALSE) > 0) %>% 
    as.data.frame() %>% 
    dplyr::as_tibble(rownames = "pop") 
  colnames(NumMLGs) <- c("pop","Num.mlgs")
  
  NumSamps <-as.data.frame(table(pop(gid))) 
  colnames(NumSamps) <- c("pop","Num.Samps")
  
  #Then pull all the pgen and psex numbers into one data.frame
  pgen.psex.by.mlg <- full_join(pgen.both.exp,
                                psex.second.both, by = "sample") %>%
    full_join(psex.multiple.both, by = "sample") %>% 
    dplyr::mutate(mlg = mll(gid), pop = pop(gid)) %>% 
    left_join(NumMLGs, by = "pop") %>% 
    left_join(NumSamps, by = "pop") %>% 
    add_count(mlg) %>% 
    arrange(n, mlg, desc(Psex.multiple))
  
}



  