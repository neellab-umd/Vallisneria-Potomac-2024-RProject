## This function runs poppr on our genclone files and then manipulates the output to create the variables we want to present.  In so doing, we delete a number of variables calculated by poppr.  Those could be added back in if desired.

## The basic poppr summary table only rarefies the estimate of eMLG.  Use this code to get rarefied values for the other statistics you need by specifying the number of samples to rarefy by the minimum sample size from each population in the n.rare argument.

## The function also cleans things up.  E.g., poppr has an unfortunate output format that puts the upper and lower CI values in one column with parentheses and a comma. So you need to get those numbers into their own columns so you can work with them. There are ensuing issues with numbers being character vectors that we have to convert.

## When you do rarefaction, poppr puts the population column in as a rowname, the function gets that back into a column

## The function also calculates the normalized effective Diversity statistics by converting Shannons (H) and Simpsons (lambda) and their CI'S into effective MLGs (i.e. exp(H) and 1/(1-lambda)

create.poppr.tables <- function(gc, group, rarefy, n.rare){

#  Set the Pop to the user specified group.

  setPop(gc)<-as.formula(paste0("~",noquote(group)))

#Run poppr
  poppr.table <- poppr(gc, minsamp = n.rare, plot = FALSE, hist = FALSE)
  poppr.table.clean <- poppr.table %>% 
    dplyr::select(Pop, N, MLG, eMLG, H, lambda) %>% 
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
    dplyr::select(-H, -lambda) %>% 
    relocate(GD, .before = GD.SW)

# Get labels from the chosen strata for each chosen group from the genclone object to add to the output
strata.labels<-as.data.frame(gc@strata) %>% 
  dplyr::select(all_of(group)) %>% 
  distinct(., .keep_all = TRUE) %>% 
  remove_rownames()

#Join the poppr data with the grouping labels.
final.poppr.table<-left_join(strata.labels,poppr.table.clean, by=structure(names = group, .Data = "Pop"))

if(rarefy == TRUE) {
   
  diversity.table<-diversity_ci(gc, 
                                rarefy = rarefy, 
                                raw = FALSE, 
                                n.rare = n.rare,
                                plot = FALSE)
  
  diversity.table.clean<-diversity.table %>% 
    tibble::rownames_to_column(var="Pop") %>% 
    dplyr::select(Pop, H.est, H.ci,lambda.est,lambda.ci) %>% 
    mutate(Eff.SW.est = exp(as.numeric(H.est)),
           Eff.Simp.est = 1/(1-as.numeric(lambda.est))) %>% 
    mutate(GD.SW.est = (Eff.SW.est-1)/(n.rare-1),
           GD.Simp.est = (Eff.Simp.est-1)/(n.rare-1)) %>% 
    mutate_all(~(gsub("\\(", "",.))) %>% #remove parentheses
    mutate_all(~(gsub("\\)", "",.))) %>% 
   #split CI columns into low and high values
    separate(H.ci, c("H.Low.CI", "H.Upper.CI"), ",")%>%
    separate(lambda.ci, c("lambda.Low.CI", "lambda.Upper.CI"), ",") %>% 
    mutate_at(vars(-Pop), ~as.numeric(.)) %>% #make all vectors except Pop numeric
    mutate(Eff.SW.Lower.CI = exp(as.numeric(H.Low.CI)),
           Eff.SW.Upper = exp(as.numeric(H.Upper.CI)),
           Eff.Simp.Lower.CI = 1/(1-as.numeric(lambda.Low.CI)),
           Eff.Simp.Upper.CI = 1/(1-as.numeric(lambda.Upper.CI))) %>% 
    dplyr::select(-starts_with("H"), -starts_with("lambda")) %>% #delete non-standardized indices
    mutate_at(vars(-Pop, -starts_with("GD")), 
              ~round(.,1)) %>% 
    mutate_at(vars(starts_with("GD")), 
              ~round(.,2)) %>% 
    dplyr::rename_at(vars(-Pop),
                     ~paste0(.x,".",n.rare)) 
 
  #Join rarefied data with the basic poppr data
  full.table<-left_join(final.poppr.table, 
                        diversity.table.clean, 
                        by=structure(names = group, .Data = "Pop"))

  return(full.table)
}
  else
  return(final.poppr.table)
}
