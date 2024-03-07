string2fasta <- function(dfWorking = dfWorking, bed_inferred_left = bed_inferred_left, bed_inferred_right = bed_inferred_right){
  
  df_deltafinal <- tibble(String = character(), 
                          fasta = character())
  
  for(i in 1:nrow(dfWorking)){  
    
    #Pull reference sequence for coordinates listed (Full Region)   #######################################Temp because we'll ultimately have the sequence...
    df_ucscOutput <- ucsc_genome_api_seq(dfWorking[[i, 3]], dfWorking[[i, 4]] - 1, dfWorking[[i, 5]])
    seq_Query <- toupper(df_ucscOutput$content$dna)
    
    #Pull reference sequence for coordinates listed (ISFG Minimum Req)
    df_ucscOutput <- ucsc_genome_api_seq(dfWorking[[i, 3]], dfWorking[[i, 6]] - 1, dfWorking[[i, 7]])
    seq_ISFG <- toupper(df_ucscOutput$content$dna)
    
    #Pull reference sequence for coordinates listed (Left Anchor)
    df_ucscOutput <- ucsc_genome_api_seq(dfWorking[[i, 3]], dfWorking[[i, 8]] - 1, dfWorking[[i, 9]])
    seq_LA <- toupper(df_ucscOutput$content$dna)
    
    #Pull reference sequence for coordinates listed (Right Anchor)
    df_ucscOutput <- ucsc_genome_api_seq(dfWorking[[i, 3]], dfWorking[[i, 10]] - 1, dfWorking[[i, 11]])
    seq_RA <- toupper(df_ucscOutput$content$dna)
    
    #Pull reference sequence for coordinates listed (Left Inferred)
    if(dfWorking[[i, 8]] < dfWorking[[i, 4]]){
      
      df_ucscOutput <- ucsc_genome_api_seq(bed_inferred_left[[i, 1]], bed_inferred_left[[i, 2]], bed_inferred_left[[i, 3]])
      seq_LI <- toupper(df_ucscOutput$content$dna)
      
    }else{}
    
    #Pull reference sequence for coordinates listed (Right Inferred)
    if(dfWorking[[i, 5]] < dfWorking[[i, 10]]) {
      
      df_ucscOutput <- ucsc_genome_api_seq(bed_inferred_right[[i, 1]], bed_inferred_right[[i, 2]], bed_inferred_right[[i, 3]])
      seq_RI <- toupper(df_ucscOutput$content$dna)
      
    }else{}
    
    #Merge inferred bases, if necessary, with query sequence
    if(dfWorking[[i, 8]] < dfWorking[[i, 4]] & dfWorking[[i, 5]] < dfWorking[[i, 10]]){
      fasta_amended <- paste0(seq_LI, seq_Query, seq_RI)
    }else{
      if(dfWorking[[i, 8]] < dfWorking[[i, 4]]){
        fasta_amended <- paste0(seq_LI, seq_Query)
      }else{
        if(dfWorking[[i, 5]] < dfWorking[[i, 10]]){
          fasta_amended <- paste0(seq_Query, seq_RI)
        }else{
          fasta_amended <- seq_Query
        }
      }
    }
    
    #Join before and after of string
    df_delta <- tibble(String = dfWorking[[i, 2]], 
                       fasta = fasta_amended) %>% 
      add_row() %>% 
      add_row() %>% 
      add_row()
    
    #Append current string to running
    df_deltafinal <- bind_rows(df_deltafinal, df_delta)
    
  }
  
  sampleName <- sub(" ", "_", paste0("output_", format(Sys.time(), "%d-%b-%Y %H.%M")))
  
  write_tsv(as_tibble(df_deltafinal$fasta), 
            file = paste0("output/", sampleName, ".fasta"))
  
  df_deltafinal <- df_deltafinal %>% 
    filter(!is.na(String))
  
  return(list(sampleName, df_deltafinal))
}

fasta2hap <- function(sampleName = sampleName, df_deltafinal = df_deltafinal){
  
  df_allSeq_final <- tibble(LocusAllele = character(), 
                            AlleleLength = character(),
                            Haplotype = character(), 
                            HaplotypeCount = numeric(),
                            HaplotypeRC = numeric())
  
  ActiveFASTA <- paste0("output/", sampleName, ".fasta")
  
  for(i in 1:nrow(df_deltafinal)){
    
    k <- (4 * i) - 4
    l <- k + 5
    m <- nrow(df_deltafinal)
    n <- m * 4
    
    #Create fasta and STRait Razor command for system call
    if(i == 1){
      
      STR8RZR_Call <- paste0("sed 5,", n, "d ", ActiveFASTA, " | ./bin/str8rzr -c db/configs/", v3_Config, " -p ", numcores)
      
    }else{
      if(i == m){
        
        STR8RZR_Call <- paste0("sed 1,", k, "d ", ActiveFASTA, " | ./bin/str8rzr -c db/configs/", v3_Config, " -p ", numcores) 
        
      }else{
        
        STR8RZR_Call <- paste0("sed ", l, ",", n, "d ", ActiveFASTA, " | sed ", 1, "," , k, "d | ./bin/str8rzr -c db/configs/", v3_Config, " -p ", numcores)
        
      }
    }
    
    #Run STRait Razor on string
    if(identical(system(STR8RZR_Call, intern = TRUE), character(0))){
      df_allSeq <- tibble(LocusAllele = "-", 
                          AlleleLength = "-",
                          Haplotype = "-", 
                          HaplotypeCount = 0,
                          HaplotypeRC = 0)
    }else{
      
      df_allSeq <- read.table(text = system(STR8RZR_Call, 
                                            intern = TRUE), 
                              sep = "\t")
    }
    
    
    names(df_allSeq) <- c("LocusAllele", "AlleleLength", "Haplotype", "HaplotypeCount", "HaplotypeRC")
    
    #Append current string to running
    df_allSeq_final <- bind_rows(df_allSeq_final, df_allSeq)
    
  }
  
  df_allSeq_final <- df_allSeq_final %>% 
    as_tibble(rownames = "Counter")
  
  df_final <- left_join(as_tibble(select(dfWorking, -String), rownames = "Counter"), df_deltafinal, by = "Counter") 
  df_final <- left_join(df_final, select(df_allSeq_final, Counter, Haplotype), by = "Counter") %>%
    select(-Counter)
  
  write_tsv(df_final, paste0("output/", sampleName, "_final.tsv"))
  
  return(df_final)
  
}