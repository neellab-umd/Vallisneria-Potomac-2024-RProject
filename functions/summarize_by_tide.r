## This function gives you the mean values for chosen population genetic statistics for each tidal regime and overall across both tidal and nontidal sites. Because it will use all numeric variables, you select the variables that you want to pass to the function. 

#M. Neel 2024.

summarize_by_tide <- function(df) {
  df %>% 
    bind_rows(., df %>% mutate(Tide = as.character("All"))) %>%
    group_by(Tide) %>% 
    dplyr::summarize(n = n(),
                     across(where(is.numeric), 
                            list(mean = ~mean(.x,na.rm = TRUE), 
                                 sd = ~sd(.x,na.rm = TRUE),
                                 min = ~min(.x,na.rm = TRUE), 
                                 max = ~max(.x,na.rm = TRUE))
                            ) 
                     ) %>% 
    as.data.frame(.) %>% 
    gather(stat, val, -c(Tide, n)) %>%
    separate(stat, into = c("variable", "stat"), sep = "_") %>%
    spread(stat, val) %>%
    dplyr::select(Tide,n,variable, min, max, mean, sd) %>% 
    print.data.frame(.,digits = 3) # reorder columns
}
