overall.poppr.table <- function(gc, group, rarefy, n.rare){

#This function runs poppr on our genclone files and then manipulates the output to create the variables we want to present only for the Total Row in the output table. So it only gives the overall pooled values, not the values for each strat In so doing, we delete a number of variables calculated by poppr.  Those could be added back in if desired.

## The basic poppr summary table only rarefies the estimate of eMLG.  Use this code to get rarefied values for the other statistics you need by specifying the number of samples to rarefy by the minimum sample size from each population in the n.rare argument.

## The function also cleans things up.  E.g., poppr has an unfortunate output format that puts the upper and lower CI values in one column with parentheses and a comma. So you need to get those numbers into their own columns so you can work with them. There are ensuing issues with numbers being character vectors that we have to convert.

## When you do rarefaction, poppr puts the population column in as a rowname, the function gets that back into a column

## The function also calculates the normalized effective Diversity statistics by converting Shannons (H) and Simpsons (lambda) and their CI'S into effective MLGs (i.e. exp(H) and 1/(1-lambda)

#  Set the Pop to the user specified group, i.e. which strata do you want to consider to be your population.

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


  return(poppr.table.clean)
}
 