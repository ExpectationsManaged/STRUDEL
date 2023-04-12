  ###########Define user agent##########
  #Identifies user to api
  ua <- user_agent("https://github.com/ExpectationsManaged/STRUDEL")

ucsc_genome_api_var <- function(chr, start, stop, 
                                  track = c("All_SNPs_(151)", "Common_SNPs_(151)", "gnomAD_v3.1")) {
  
    tracks <- tibble(shortLabel = c("snp151", "snp151Common", "gnomadGenomesVariantsV3_1"), longLabel = c("All_SNPs_(151)", "Common_SNPs_(151)", "gnomAD_v3.1"))
  
    activeTrack <- tracks %>%
      filter(longLabel == track)
  
    activeTrack <- activeTrack[[1, 1]]
    
    path <- paste0("/getData/track?genome=hg38;track=", activeTrack, ";chrom=chr", str_remove({chr}, "chr"),";start=", {start}, ";end=", {stop})
    url <- paste0("https://api.genome.ucsc.edu", path = path)
  
    #Retrieve url
    resp <- GET(url, ua)
    
    #Check type and confirm JSON
    if (http_type(resp) != "application/json") {
      stop("API did not return json", call. = FALSE)
    }
    
    #Parse JSON result from GET url call
    parsed <- jsonlite::fromJSON(content(resp, "text"), simplifyVector = FALSE)
    
    #Return error information from GET url call
    if (http_error(resp)) {
      stop(
        sprintf(
          "GitHub API request failed [%s]\n%s\n<%s>", 
          status_code(resp),
          parsed$message,
          parsed$documentation_url
        ),
        call. = FALSE
      )
    }
    
    #Create parsed S3 object
    structure(
      list(
        content = parsed,
        path = path,
        response = resp
      ),
      class = "ucsc_genome_api_var"
    )
  }