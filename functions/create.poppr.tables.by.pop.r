create.poppr.tables.by.pop <- function(gc,rarefy, n.rare){

  #This function runs poppr on our genclone files and then manipulates the output to create the variables we want to present.  In so doing, we delete a number of variables calculated by poppr.  Those could be added back in if desired.
#
## The basic poppr summary table only rarefies the estimate of eMLG.  Use this code to get rarefied values for the other statistics you need by specifying the number of samples to rarefy by the minimum sample size from each population.
#
## The function also cleans things up.  E.g., poppr has an unfortunate output format that puts the upper and lower CI values in one column with parentheses and a comma. So you need to get those numbers into their own columns so you can work with them. There are ensuing issues with numbers being character vectors that we have to convert.
## When you do rarefaction, poppr puts the population column in as a rowname, the function gets that back into a column
  
## The function also calculates the normalized effective Diversity statistics by converting Shannons (H) and Simpsons (lambda) and their CI'S into effective MLGs (i.e. exp(H) and 1/(1-lambda)

  poppr.table <- poppr(gc, minsamp = n.rare, plot = FALSE, hist = FALSE)
  poppr.table.clean <- poppr.table %>% 
    select(Pop, N, MLG, eMLG, H, lambda) %>% 
    mutate(GD = (MLG-1)/(N-1)) %>% 
    mutate(Eff.SW = exp(as.numeric(H)),
           Eff.Simp = 1/(1-as.numeric(lambda))) %>% 
    mutate(GD.SW = (Eff.SW-1)/(N-1),
           GD.Simp = (Eff.Simp-1)/(N-1)) %>% 
    mutate_at(vars(-Pop), ~as.numeric(.)) %>% 
    mutate_at(vars(-Pop, -N, -MLG, -starts_with("GD")), 
              ~round(.,1)) %>% 
    mutate_at(vars(starts_with("GD")), 
              ~round(.,2)) %>% 
    select(-H, -lambda) %>% 
    relocate(GD, .before = GD.SW)

# Get strata labels for each population to add to the output
strata.labels<-as.data.frame(gc@strata) %>% 
  distinct(NewPop, .keep_all = TRUE) %>% 
  remove_rownames()

final.poppr.table<-left_join(strata.labels,poppr.table.clean, by=c("NewPop" = "Pop"))

if(rarefy == TRUE) {
   
  diversity.table <- diversity_ci(gc, 
                                rarefy = rarefy, 
                                raw = FALSE, 
                                n.rare = n.rare,
                                plot = FALSE)
  
  diversity.table.clean <- diversity.table %>% 
    tibble::rownames_to_column(var="Pop") %>% 
    select(Pop, H.est, H.ci,lambda.est,lambda.ci) %>% 
    mutate(Eff.SW.est = exp(as.numeric(H.est)),
           Eff.Simp.est = 1/(1-as.numeric(lambda.est))) %>% 
    mutate(GD.SW.est = (Eff.SW.est-1)/(n.rare-1),
           GD.Simp.est = (Eff.Simp.est-1)/(n.rare-1)) %>% 
    mutate_all(~(gsub("\\(", "",.))) %>% 
    mutate_all(~(gsub("\\)", "",.))) %>% 
    #split CI columns into a low and high value
    separate(H.ci, c("H.Low.CI", "H.Upper.CI"), ",")%>%
    separate(lambda.ci, c("lambda.Low.CI", "lambda.Upper.CI"), ",") %>% 
    mutate_at(vars(-Pop), ~as.numeric(.)) %>%
    mutate(Eff.SW.Lower.CI = exp(as.numeric(H.Low.CI)),
           Eff.SW.Upper = exp(as.numeric(H.Upper.CI)),
           Eff.Simp.Lower.CI = 1/(1-as.numeric(lambda.Low.CI)),
           Eff.Simp.Upper.CI = 1/(1-as.numeric(lambda.Upper.CI))) %>% 
    select(-starts_with("H"), -starts_with("lambda")) %>% 
    mutate_at(vars(-Pop, -starts_with("GD")), 
              ~round(.,1)) %>% 
    mutate_at(vars(starts_with("GD")), 
              ~round(.,2)) %>% 
    dplyr::rename_at(vars(-Pop),
                     ~paste0(.x,".",n.rare)) 
 
  full.table<-left_join(final.poppr.table, 
                        diversity.table.clean, 
                        by=c("NewPop" = "Pop"))

  return(full.table)
}
  else
  return(final.poppr.table)
}
