### Hypergeometric probability of finding overlapping features in two random samples
### from a single universe of features

# set_a {character} Vector of features identified in set A
# set_b {character} Vector of features identified in set B
# universe {character} Vector of features identified in universe

sig_overlap <- function(set_a, set_b, universe){
  if(!all(set_a %in% universe)) stop('Not all features in set A is found in the universe of features.')
  if(!all(set_b %in% universe)) stop('Not all features in set B is found in the universe of features.')
  
  count_univ <- length(universe)
  count_a <- length(set_a)
  count_b <- length(set_b)
  overlaps <- intersect(set_a, set_b)
  count_overlap <- length(overlaps)
  
  p_val <- 1 - phyper(count_overlap, count_b, count_univ - count_b, count_a)
  return(p_val)
}