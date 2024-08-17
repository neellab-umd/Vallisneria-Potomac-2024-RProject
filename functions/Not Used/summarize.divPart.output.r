#This function uses another custom function summarize.divPart.metric.pop.mean() to calculate the mean value of each metric in a list of metrics.  It then combines those mean values from each list and puts them into a useful dataframe that yields a table with mean values for each site that we can combine with the within population diversity statistics.

#First source the custom function partition.metric.pop.mean() to calculate the mean of pairwise values of each metric for each population. fastDivpart outputs Metrics outputs in a list in which each statistic is an element of the list.

  summarize.divPart.output <- function(partition.list) {

#Source the custom function
if(!exists("partition.metric.pop.mean", mode="function")) source(here::here("functions", "summarize.divPart.metric.pop.mean.r"))
    
  divPart.means <-
    purrr::map(partition.list$pairwise,
               partition.metric.pop.mean) %>% 
    #add column with the metric name you get from the list name
    dplyr::bind_rows(.id = "metric.name") %>% 
    #convert to wide data with metric in column with the metric name
    tidyr::pivot_wider(names_from = "metric.name", values_from = "metric") %>% 
    #get rid of the underscore and number genepop format uses in the population name
    tidyr::separate(Pop, c("Pop"), sep = "_", extra='drop') %>% 
    dplyr::mutate_if(is.numeric, round, 3)
  
    return(divPart.means)
  
  }
  
  
  