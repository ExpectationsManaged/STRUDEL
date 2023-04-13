getBED <- function(path2File){
  
  ###########Load bed file##########
  
  auSTRs_bed <- read_delim(path2File, 
                       delim = "\t", 
                       escape_double = FALSE, 
                       col_names = FALSE, 
                       trim_ws = TRUE)|> 
    dplyr::rename(Chr = X1, 
                  Start = X2, 
                  Stop = X3, 
                  Locus = X4)
  
  ###########Compile list of loci##########
  
  locus_list <- auSTRs_bed %>%
    distinct(Locus)
  
  return(list(auSTRs_bed, locus_list))
}
