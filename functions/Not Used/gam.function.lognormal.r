gam.function.lognormal <- function( data, yvar, xvar,) {
  gamlss(yvar ~ xvar, 
           sigma.formula = ~xvar,
           data = data, 
           gamlss.family = LNO(mu.link = "log", sigma.link = "log")) %>%  
    broom.mixed::tidy(.x)                                                                                  }

#assign(paste0("gam.sum", yvar,".by.", xvar))
#write_csv(gam.sum.lognorm.mean.NumAllStems.Within.Site.before.after, file = here("output.data", paste0("gam.sum.lognorm.mean.NumAllStems.Within.Site.before.after.csv")))
