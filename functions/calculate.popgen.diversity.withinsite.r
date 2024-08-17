## This function calculates all the basic within-site population genetic diversity statistics using adegnet, pegas, and hierfstat.

## M. Neel June 2024.

calculate.popgen.diversity.withinsite <- function(genind) {
  
  sample.size<-summary(genind)$n.by.pop
  
  Alleles<-summary(genind)$pop.n.all
  
  ## The minimum number of alleles is set manually to the number at Pennyfield Lock.  
  
  Richness <- hierfstat::allelic.richness(hierfstat::genind2hierfstat(genind))
  
  loc.n.all.by.pop <- sapply(seppop(genind, drop=TRUE), function(pop) pop@loc.n.all) 
  PercPolymorphic <- colSums(loc.n.all.by.pop > 1)/10  
  
  Hobs <- t(sapply(seppop(genind), function(ls) summary(ls)$Hobs))
  Hobs.pop <- apply(Hobs, MARGIN = 1, FUN = mean)
  
  Hexp <- t(sapply(seppop(genind), function(ls) summary(ls)$Hexp))
  Hexp.pop <- apply(Hexp, MARGIN = 1, FUN = mean) 
  
  # Function to calculate mean FIS for a single population
  calculate_mean_FIS <- function(genind) {
    pop_hier <- hierfstat::genind2hierfstat(genind)
    FIS_values <- hierfstat::basic.stats(pop_hier)$Fis
    mean(FIS_values, na.rm = TRUE)
  }
  
  # Function to perform bootstrapping and extract CI limits
  bootstrap_FIS_CI <- function(genind) {
    pop_hier <- hierfstat::genind2hierfstat(genind)
    boot_results <- hierfstat::boot.ppfis(pop_hier, nboot = 1000)
    # Extract the lower and upper confidence intervals
    fis_ci <- boot_results$fis.ci
    return(fis_ci)
  }
  
  # Calculate mean FIS for each population
  mean_FIS_pop <- sapply(seppop(genind), calculate_mean_FIS)
  
  # Perform bootstrap estimates and get CI for each population
  boot_FIS_CI_list <- lapply(seppop(genind), bootstrap_FIS_CI)
  
  # Extract lower and upper CI values
  lower_CI <- sapply(boot_FIS_CI_list, function(ci) ci$ll[1])
  upper_CI <- sapply(boot_FIS_CI_list, function(ci) ci$hl[1])
  
  # Combine the results into a data frame
  within.pop.popgen.diversity.summary <- data.frame(Pop = names(Hobs.pop),
                                                 N.Potomac.No.Reps = sample.size,
                                                 PercPolymorphic = PercPolymorphic,
                                                 A.All = as.numeric(Alleles),
                                                 Ap = Alleles/10,
                                                 Ar = round(colMeans(Richness$Ar), 2),
                                                 Ho = round(Hobs.pop,2),
                                                 He = round(Hexp.pop,2),
                                                 Fis = round(mean_FIS_pop,2),
                                                 Lower_Fis_CI = round(lower_CI,2),
                                                 Upper_Fis_CI = round(upper_CI,2)) %>% 
                                          dplyr::mutate(Fis = na_if(Fis,NaN),
                                                 Lower_Fis_CI = na_if(Lower_Fis_CI,-Inf),
                                                 Upper_Fis_CI = na_if(Upper_Fis_CI,-Inf)) 
                                                   

  return(within.pop.popgen.diversity.summary)
  
  
}

