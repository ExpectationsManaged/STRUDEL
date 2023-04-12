###########Biblioteca##########

{
library(rsnps)
library(tidyverse)
library(httr)     
library(jsonlite)
library(reticulate)
}



###########Load functions###########
source(getBED.r)
source(returnSeq.r)
source(returnVar.r)



###########Load bed file###########
list_bed <- getBED(file.choose())

	bed_loci <- list_bed[1]
	list_loci <- list_bed[2]



##########Get sequences from UCSC api using bed file##########
parsed_UCSC_JSON_bed_seq <- enframe(unlist(pmap(bed_loci[1:3], 
                                                ~returnSeq({..1}, {..2}, {..3})))) |> 
  pivot_wider(names_from = name, 
              values_from = value, 
              values_fn = list(value = list)) |>
  select(content.chrom, 
         content.start, 
         content.end, 
         content.dna) |> 
  unchop(everything())

##########Capitalize on things##########
parsed_UCSC_JSON_bed_seq$content.dna <- toupper(parsed_UCSC_JSON_bed_seq$content.dna)



###########Parse JSON file for lots of loci##########
parsed_UCSC_JSON_bed_var <- pmap(bed_loci[1:3], 
                                 ~returnVar({..1}, {..2}, {..3}, 
                                                      track = "gnomAD_v3.1"))

###########Create data frame of parsed JSON##########
contents <- lapply(parsed_UCSC_JSON_bed_var, function(x) data.frame(x$content))

###########Bind data by row and fill missing data##########
parsed_UCSC_JSON_bed_var <- do.call(plyr::rbind.fill, contents) |>  
  select(chrom, start, end, 
         contains("Start"), 
         contains("End"), 
         contains("ref"), 
         contains("rsId"), 
         contains("alt"), 
         contains("AF"), 
         contains("AC"), 
         contains("strand")) |>
  select(-contains("startPos"))

###########Replace periods with underscores because pivot longer hates me  ##########
names(parsed_UCSC_JSON_bed_var) <- gsub(".", "_", 
                                        names(parsed_UCSC_JSON_bed_var), 
                                        fixed = TRUE)

###########Replace underscores  with periods because spite  ##########
names(parsed_UCSC_JSON_bed_var) <- gsub("3_1", "3.1", 
                                               names(parsed_UCSC_JSON_bed_var), 
                                               fixed = TRUE)

###########Pivot table to long form and, optionally, filter out extremely low-freq variants##########
parsed_UCSC_JSON_bed_var <- parsed_UCSC_JSON_bed_var |>  
  pivot_longer(cols = contains("_"), 
               names_to = c(".value", "set"), 
               names_prefix = "gnomadGenomesVariantsV3.1_", 
               names_sep = "_", 
               values_drop_na = FALSE) #|>  
  # filter(avHet >= 0.001)

###########Replace first observation artifact ##########
parsed_UCSC_JSON_bed_var$set <- replace_na(parsed_UCSC_JSON_bed_var$set, 0)
