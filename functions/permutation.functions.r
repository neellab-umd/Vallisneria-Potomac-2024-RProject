## Functions for permutation testing for differences among tidal regimes.

## M. Neel 2024

#Function for calculating differences by Tidal Regime. Used for the observed data and all the permutations.

calculate_diffs <- function(perm) {
  df <- as.data.frame(perm)
  df %>%
    group_by(Tide) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")  %>% 
    pivot_longer(cols=c(-Tide),names_to="Original_Vars") %>%
    pivot_wider(names_from=c(Tide)) %>% 
    mutate(Difference = abs(Tidal - NonTidal)) %>% 
    dplyr::select(-c(Tidal, NonTidal)) %>% 
    pivot_wider(names_from = Original_Vars, values_from = Difference) 
}


# Function that calculates the p-value for a single variable

calculate_p_value <- function(observed_var,
                              observed_val,
                              simulated_values) {
  num_extreme_values <- sum(simulated_values >= observed_val)
  p_value <- num_extreme_values / length(simulated_values)
  return(tibble(variable = observed_var, p_value = p_value))
}