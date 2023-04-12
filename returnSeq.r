
###########Identify user to api##########
ua <- user_agent("https://github.com/ExpectationsManaged/STRUDEL")

  ###########Return sequence string##########
  ucsc_genome_api_seq <- function(chr, start, stop) {
    ############
    # url <- modify_url("https://api.genome.ucsc.edu", path = path)
    
    # chr <- str_remove(bed[[1,1]], "chr")
    # start <- bed[[1,2]]
    # stop <- bed[[1,3]]
    # browser()
    #Construct url for api call
    
    #################
    path <- paste0("/getData/sequence?genome=hg38;chrom=chr", 
                   str_remove({chr}, "chr"),
                   ";start=", {start}, 
                   ";end=", {stop})
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
      class = "ucsc_genome_api_seq"
    )
  }