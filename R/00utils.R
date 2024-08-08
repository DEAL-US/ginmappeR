utils::globalVariables('cardPath')

.checkNotNull <- function(variable, message="The given ID is "){
    if(identical(variable, "")){
        .terminateProgressBar()
        stop(message, "empty")
    }else if(is.null(variable)){
        .terminateProgressBar()
        stop(message, "NULL")
    }else if(is.nan(variable)){
        .terminateProgressBar()
        stop(message, "NaN")
    }else if(is.na(variable)){
        .terminateProgressBar()
        stop(message, "NA")
    }
}

.checkNCBIProteinIdExists <- function(proteinId){
    .checkNotNull(proteinId, 'The given NCBI Protein ID is ')
    # Example id: ADH03009.1
    # [ACCN] means we look for an Accesion Number OR [UID] Unique number assigned to publication
    invisible(capture.output(query <- entrez_search(db='protein', term=paste(proteinId,'[ACCN] OR ', proteinId, '[UID]',sep=''))))
    if(!length(query[['ids']])>0) {
        .terminateProgressBar()
        stop('The given ID "', proteinId, '" is not registered in NCBI Protein database')
    }
}

.checkNCBINucleotideIdExists <- function(nucleotideId){
    .checkNotNull(nucleotideId, 'The given NCBI Nucleotide ID is ')
    # Example id: HM036080.1
    invisible(capture.output(query <- entrez_search(db='nucleotide', term=paste(nucleotideId,'[ACCN] OR ', nucleotideId, '[UID]',sep=''))))
    if(!length(query[['ids']])>0) {
        .terminateProgressBar()
        stop('The given ID "', nucleotideId, '" is not registered in NCBI Nucleotide database')
    }
}

.checkNCBIGeneIdExists <- function(geneId){
    .checkNotNull(geneId, 'The given NCBI Gene ID is ')
    # Example id: WP_001082319
    invisible(capture.output(query <- entrez_search(db='gene', term=paste(geneId,'[ACCN] OR ', geneId, '[UID]',sep=''))))
    if(!length(query[['ids']])>0) {
        .terminateProgressBar()
        stop('The given ID "', geneId, '" is not registered in NCBI Gene database')
    }
}

.checkUniProtIdExists <- function(upId){
    .checkNotNull(upId, 'The given UniProt ID is ')
    # Example id not registered: P0ZTH2
    # Example id not valid: P0ZTUH2
    # Example valid id: P0DTH6
    tryCatch(
        {
            invisible(capture.output(queryUniProt(
                query = paste('accession:', upId, sep=''),
                fields = c('accession'),
                collapse = ' OR '
            )))
        },
        error = function(e) {
            .terminateProgressBar()
            stop('The given ID "', upId, '" is not registered in UniProt database')
        })
}


.checkKEGGIdExists <- function(keggId){
    .checkNotNull(keggId, 'The given KEGG ID is ')

    # Search ID in KEGG genes database
    # query <- keggFind('genes', keggId)
    #
    # # Check if there is a coincidence and if it corresponds to an unique gene
    # if(!length(query)==1){
    # .terminateProgressBar()
    # stop('The given ID "', keggId, '" is not registered in KEGG genes database')
    # }
}

.checkCARDIdExists <- function(cardId){
    .checkNotNull(cardId, 'The given CARD ARO accession is ')

    checkedCardId <- if (grepl("^ARO:", cardId)) cardId else paste0("ARO:", cardId)

    aroIndex <- .loadAROIndex()

    if(!checkedCardId %in% aroIndex$ARO.Accession) {
        .terminateProgressBar()
        stop('The given ID "', cardId, '" is not registered in CARD database')
    }
}

.checkClusterIdentity <- function(clusterIdentity){
    if(!clusterIdentity %in% c('1.0','0.9','0.5')){
        .terminateProgressBar()
        stop('The requested cluster identity "', clusterIdentity, '" is not valid')
    }
}

.checkBoolean <- function(value, name){
    if(!value %in% c(TRUE, FALSE)){
        .terminateProgressBar()
        stop(name, ' variable value must be TRUE or FALSE, value given: ', value)
    }
}

.mergeNamedLists <- function(list1,list2){
    for (key in names(list2)) {
        if (key %in% names(list1)) {
            list1[[key]] <- c(list1[[key]], list2[[key]])
        } else {
            list1[[key]] <- list2[[key]]
        }
    }
    return(list1)
}

.checkIdInNCBIDataBase <- function(ncbiId, ncbiDB) {
    result<- tryCatch(
        switch (ncbiDB,
                'protein' = .checkNCBIProteinIdExists(ncbiId),
                'nucleotide' = .checkNCBINucleotideIdExists(ncbiId),
                'gene' = .checkNCBIGeneIdExists(ncbiId),
        ),
        error = function(e){return(FALSE)}
    )
    if(is.null(result)){result <- TRUE}
    return(result)
}

.testEquals <- function(testFunction, result){
    tryCatch(
        {RUnit::checkEquals(testFunction, result)},
        error = function(e){
            message('API(s) connection failed:')
            message(e)
        }
    )
}

.suppressOutput <- function(expr) {
    capture.output(suppressWarnings(suppressMessages({
        result <- expr
    })), file = NULL)
    return(result)
}

#####################################
# CARD database retrieval functions #
#####################################

.downloadAndExtractCARD <- function(){
    if(is.null(getOption("cardPath"))){
        options("cardPath" = tempdir())
    }

    localRelativeFilenames <- c(paste('card-data.tar.bz2-', format(Sys.time(), "%a %m-%d-%Y %H-%M-%S"),sep=''))
    url <- c('https://card.mcmaster.ca/latest/data')

    # Delete previous versions
    message('Deleting previous versions of CARD if any.', "\n")
    unlink(paste(getOption("cardPath"),'/card-data*',sep=''), recursive = TRUE)

    # Download CARD zip
    message('Downloading latest version', "\n")
    download.file(url, paste(getOption("cardPath"),'/',localRelativeFilenames,sep=''), quiet = TRUE, mode = 'wb')
    # Extract it
    message('Extracting database', "\n")
    zipF<- paste(getOption("cardPath"),'/',localRelativeFilenames, sep='')
    # Unzip the compressed file
    untar(zipF, exdir = paste(getOption("cardPath"),'/card-data',sep=''))
}

updateCARDDataBase <- function(){
    tryCatch(
        {
            message('Updating CARD database data...')
            .downloadAndExtractCARD()
            message('CARD database downloaded successfully!\nLocated at ', getOption("cardPath"),'/card-data')
            .getCARDVersionAux()
        },
        error = function(e) {return(message('The download of CARD database failed, please try again.\nIf the error persists try using changeCARDPath() to refresh CARD folder location.\nIf not, consider downloading and extracting it manually. Then, use changeCARDPath() function to indicate the folder where it is located.', "\n"))},
        warning = function(e) {return(message('The download of CARD database failed, please try again.\nIf the error persists try using changeCARDPath() to refresh CARD folder location.\nIf not, consider downloading and extracting it manually. Then, use changeCARDPath() function to indicate the folder where it is located.', "\n"))}
    )
}

.checkIfCARDIsDownloaded <- function(){
    if(is.null(getOption("cardPath"))){
        updateCARDDataBase()
    }else{
        if(!file.exists(paste(getOption("cardPath"),'/card-data/aro_index.tsv',sep=''))){
            updateCARDDataBase()
        }else{
            # Get date and time
            downloadDate <- strsplit(list.files(getOption("cardPath"), pattern = 'card-data.tar.bz2-')[[1]],'bz2-')[[1]][[2]]
            # Format it correctly
            downloadDate <- format(strptime(downloadDate, format = "%a %m-%d-%Y %H-%M-%S"), format = "%a %m/%d/%Y %H:%M:%S")
            message(sprintf('Using a CARD database version downloaded on %s, please consider updating it with updateCARDDataBase() function.',
                            downloadDate))
        }
    }
}

changeCARDPath <- function(path=tempdir()){
    if(dir.exists(path)){
        if(!is.null(options("cardPath")[[1]])){
            unlink(paste(options("cardPath")[[1]],'/card-data*',sep=''), recursive = TRUE)
        }
        options("cardPath" = path)
        updateCARDDataBase()
    }else{
        warning('Path ',path,' does not exist')
    }
}

.loadAROIndex <- function(){
    return(read.csv(paste(options("cardPath")[[1]],'/card-data/aro_index.tsv',sep=''), sep='\t'))
}

# This function translates a bunch of NCBI ids to CARD
.auxNCBI2CARD <- function(translations, ncbiDB, aroIndex){
    if(identical(ncbiDB, 'protein')){
        translations <- lapply(translations, FUN = function(x) unlist(sapply(x, USE.NAMES = FALSE, FUN = function(auxId){
            auxId <- strsplit(auxId,'.', fixed = TRUE)[[1]][[1]]
            return(aroIndex[gsub("\\..*", "", aroIndex$Protein.Accession) == auxId, "ARO.Accession"])
        })))
    }else{
        translations <- lapply(translations, FUN = function(x) unlist(sapply(x, USE.NAMES = FALSE, FUN = function(auxId){
            auxId <- strsplit(auxId,'.', fixed = TRUE)[[1]][[1]]
            return(aroIndex[gsub("\\..*", "", aroIndex$DNA.Accession) == auxId, "ARO.Accession"])
        })))
    }
    return(translations)
}

###########################################
#  Functions for caching and vectorizing  #
###########################################

# Function for cache management when errors or warnings do happen
.cacheErrorHandler <- function(res, raiseError = FALSE){
    if(is.list(res) && !is.data.frame(res) && is.null(names(res)) && length(res)>0 && is.na(res[[1]])){
        # If an error was produced, manually delete wrong cache values by comparing the keys
        # before the function call and after it
        for (x in setdiff(.getCache()$keys(), res[[2]])){
            .getCache()$remove(x)
        }
        if(raiseError){
            stop('Error ocurred')
        }else{
            return(NULL)
        }
    }else{
        # If no error, return the result
        return(res)
    }
}

.retryHandler <- function(f, ...){
    res <- NULL
    counter <- 0
    while (is.null(res) && counter<20) {
        if(counter == 19){
            res <- .cacheErrorHandler(f(...))
        }else{
            res <- suppressMessages(.cacheErrorHandler(f(...)))
        }
        counter <- counter + 1
        Sys.sleep(0.1)
    }
    return(res)
}

# .resultParser <- function(result){
#     if (length(result) == 1) {
#         return(result[[1]])
#     } else {
#         return(result)
#     }
# }

.replaceCharacter0WithNULL <- function(lst) {
    if (is.list(lst)) {
        lapply(lst, .replaceCharacter0WithNULL)
    } else if (identical(lst, character(0))) {
        NULL
    } else {
        lst
    }
}

.resultParser <- function(inputVector, exhaustiveMapping, result, isDataFrame = FALSE){
    # Unwrap result from Vectorize
    if(length(result)==1){
        result <- result[[1]]
    }

    # Check if input is a vector of IDs and if it has several ids
    # If it is character(0)
    if ((is.vector(inputVector) && !is.list(inputVector) && length(inputVector)==0)){
        if(isTRUE(exhaustiveMapping)){result <- list()
        }else{result <- character(0)}
    }
    # If it has more than 1 value
    else if (is.vector(inputVector) && !is.list(inputVector) && length(inputVector)>1) {

        # For simple mapping functions, just merge results in a vector
        if(is.null(exhaustiveMapping)){
            result[sapply(result, is.null)] <- NA
            result <- unlist(result, use.names = FALSE)
        }else{
            # For functions that may obtain several ids for mapping one id

            # Check that the result is not a dataframe so we dont flat it
            if(!isDataFrame){
                # Depending of exhaustiveMapping, we return a list or a vector and treat missing values consequently, NA for vector, NULL for list
                if(isFALSE(exhaustiveMapping)){
                    result[sapply(result, is.null)] <- NA
                    result <- lapply(result, function(a) if (is.vector(a) && !is.data.frame(a) && is.null(names(a)) && length(a) == 0) NA else a)
                    result <- unlist(result, use.names = TRUE)
                }else{
                    result <- .replaceCharacter0WithNULL(result)
                    result <- lapply(result, function(a) if (is.vector(a) && !is.data.frame(a) && is.null(names(a)) && length(a) == 0) NULL else a)
                }
            }
        }
    }else{# If the input vector has only 1 id

        # If error or no translation, return NA if exhaustiveMapping = FALSE, list(NULL) otherwise
        if(is.null(result) || (is.vector(result) && !is.data.frame(result) && is.null(names(result)) && length(result) == 0)){
            result <- NA
            if(isTRUE(exhaustiveMapping)){
                result <- list(NULL)
            }
        }else{ # If there's no error and there are translations
            if(!isDataFrame && isTRUE(exhaustiveMapping)){
                # Get the result in a list if exhaustiveMapping is set to TRUE
                result <- .replaceCharacter0WithNULL(result)
                result <- list(result)
            }
        }
    }
    return(result)
}

.errorMessageHandler <- function(e){
    if(grepl("given ID", e)){
        message(e$message)
    }else{
        # message('One of the web APIs failed, please try again.\n')
    }
    # message('One of the web APIs failed, please try again.\n', ifelse(grepl("given ID", e), e, ""), "\n")
    return(list(NA, .getCache()$keys()))
}

.warningMessageHandler <- function(e){
    # strsplit(conditionMessage(e), '\n', fixed=T)
    # warning('One of the web APIs failed, please try again.\n', "\n", call. = FALSE, noBreaks. = TRUE, immediate. = TRUE)
    # message('One of the web APIs failed, please try again.\n')
    return(list(NA, .getCache()$keys()))
}

.getCache <- function(){
    # LRU cache of 4gb of maximum size and a day of maximum age of a translation
    return(cachem::cache_disk(paste(tempdir(),'/cache',sep=''),max_size = 4096 * 1024^2, max_age = 86400, evict='lru'))
}

.createCache <- function(){

    # --- Between NCBI databases ---
    # .getNCBIDatabasesLinks <<- memoise::memoise(.getNCBIDatabasesLinks, cache = cm)
    .getNCBIGene2NCBIProtein <<- memoise::memoise(.getNCBIGene2NCBIProtein, cache = .getCache())
    .getNCBIProtein2NCBIGene <<- memoise::memoise(.getNCBIProtein2NCBIGene, cache = .getCache())
    .getNCBIProtein2NCBINucleotide <<- memoise::memoise(.getNCBIProtein2NCBINucleotide, cache = .getCache())
    .getNCBINucleotide2NCBIProtein <<- memoise::memoise(.getNCBINucleotide2NCBIProtein, cache = .getCache())
    .getNCBIGene2NCBINucleotide <<- memoise::memoise(.getNCBIGene2NCBINucleotide, cache = .getCache())
    .getNCBINucleotide2NCBIGene <<- memoise::memoise(.getNCBINucleotide2NCBIGene, cache = .getCache())

    # --- NCBI to UniProt ---
    # .getNCBI2UniProtDT <<- memoise::memoise(.getNCBI2UniProtDT, cache = .getCache())
    .getNCBIIdenticalProteins <<- memoise::memoise(.getNCBIIdenticalProteins, cache = .getCache())
    # .getNCBI2UniProtBatch <<- memoise::memoise(.getNCBI2UniProtBatch, cache = .getCache())
    # .getNCBI2UniProtIP <<- memoise::memoise(.getNCBI2UniProtIP, cache = .getCache())
    .getNCBIProtein2UniProt <<- memoise::memoise(.getNCBIProtein2UniProt, cache = .getCache())
    .getNCBINucleotide2UniProt <<- memoise::memoise(.getNCBINucleotide2UniProt, cache = .getCache())
    .getNCBIGene2UniProt <<- memoise::memoise(.getNCBIGene2UniProt, cache = .getCache())

    # --- NCBI to KEGG ---
    # .getNCBI2KEGGDT <<- memoise::memoise(.getNCBI2KEGGDT, cache = .getCache())
    # .getNCBI2KEGGTUP <<- memoise::memoise(.getNCBI2KEGGTUP, cache = .getCache())
    # .getNCBI2KEGG <<- memoise::memoise(.getNCBI2KEGG, cache = .getCache())
    .getNCBIProtein2KEGG <<- memoise::memoise(.getNCBIProtein2KEGG, cache = .getCache())
    .getNCBINucleotide2KEGG <<- memoise::memoise(.getNCBINucleotide2KEGG, cache = .getCache())
    .getNCBIGene2KEGG <<- memoise::memoise(.getNCBIGene2KEGG, cache = .getCache())

    # --- NCBI to CARD ---
    .getNCBIProtein2CARD <<- memoise::memoise(.getNCBIProtein2CARD, cache = .getCache())
    .getNCBINucleotide2CARD <<- memoise::memoise(.getNCBINucleotide2CARD, cache = .getCache())
    .getNCBIGene2CARD <<- memoise::memoise(.getNCBIGene2CARD, cache = .getCache())

    # --- CARD functions caching ----
    .getCARD2NCBIAux <<- memoise::memoise(.getCARD2NCBIAux, cache = .getCache())
    .getCARD2NCBIGene <<- memoise::memoise(.getCARD2NCBIGene, cache = .getCache())
    .getCARD2UniProt <<- memoise::memoise(.getCARD2UniProt, cache = .getCache())
    .getCARD2KEGG <<- memoise::memoise(.getCARD2KEGG, cache = .getCache())

    # --- KEGG functions caching ----
    .getKEGG2UniProt <<- memoise::memoise(.getKEGG2UniProt, cache = .getCache())
    # .getKEGG2NCBIDT <<- memoise::memoise(.getKEGG2NCBIDT, cache = .getCache())
    # .getKEGG2NCBITUP <<- memoise::memoise(.getKEGG2NCBITUP, cache = .getCache())
    # .getKEGG2NCBI <<- memoise::memoise(.getKEGG2NCBI, cache = .getCache())
    .getKEGG2NCBIProtein <<- memoise::memoise(.getKEGG2NCBIProtein, cache = .getCache())
    .getKEGG2NCBINucleotide <<- memoise::memoise(.getKEGG2NCBINucleotide, cache = .getCache())
    .getKEGG2NCBIGene <<- memoise::memoise(.getKEGG2NCBIGene, cache = .getCache())
    # .getKEGG2CARDexhaustiveMappingAux <<- memoise::memoise(.getKEGG2CARDexhaustiveMappingAux, cache = .getCache())
    .getKEGG2CARD <<- memoise::memoise(.getKEGG2CARD, cache = .getCache())

    # --- UniProt to KEGG ----
    # .getUniProt2KEGGDT <<- memoise::memoise(.getUniProt2KEGGDT, cache = .getCache())
    .getUniProtSimilarGenes <<- memoise::memoise(.getUniProtSimilarGenes, cache = .getCache())
    # .getUniProt2KEGGSGT <<- memoise::memoise(.getUniProt2KEGGSGT, cache = .getCache())
    .getUniProt2KEGG <<- memoise::memoise(.getUniProt2KEGG, cache = .getCache())

    # --- UniProt to NCBI ----
    # .getUniProt2NCBIDT <<- memoise::memoise(.getUniProt2NCBIDT, cache = .getCache())
    # .getUniProt2NCBISGT <<- memoise::memoise(.getUniProt2NCBISGT, cache = .getCache())
    # .getUniProt2NCBI <<- memoise::memoise(.getUniProt2NCBI, cache = .getCache())
    .getUniProt2NCBIProtein <<- memoise::memoise(.getUniProt2NCBIProtein, cache = .getCache())
    .getUniProt2NCBINucleotide <<- memoise::memoise(.getUniProt2NCBINucleotide, cache = .getCache())
    .getUniProt2NCBIGene <<- memoise::memoise(.getUniProt2NCBIGene, cache = .getCache())

    # --- UniProt to CARD ----
    # .getUniProt2CARDexhaustiveMappingAux <<- memoise::memoise(.getUniProt2CARDexhaustiveMappingAux, cache = .getCache())
    .getUniProt2CARD <<- memoise::memoise(.getUniProt2CARD, cache = .getCache())
}

######################################
# Progress bar management functions  #
######################################

# Function to initialize the progress bar
.initializeProgressBar <- function(total_steps) {
    .global_step_counter <<- 0
}

# Function to clear the console line after updating the progress bar
.clearConsoleLine <- function() {
    # Redirect output to the console
    sink(file = stdout(), type = "output")
    sink(file = stderr(), type = "message")
    cat("\r", strrep(" ", getOption("width")), "\r")
    flush.console()
        # Reset the output redirection
    if (sink.number(type = "output") > 0) {
        sink(type = "output")
    }
    if (sink.number(type = "message") > 0) {
        sink(type = "message")
    }
}

# Function to update the progress bar
.updateProgressBar <- function(message) {
    .global_step_counter <<- .global_step_counter + 1
    .clearConsoleLine()

    # Calculate the number of steps left
    fixed_size <- 30
    current_step <- .global_step_counter %% fixed_size
    stepsLeft <- fixed_size - current_step

    # Redirect output to the console
    sink(file = stdout(), type = "output")
    sink(file = stderr(), type = "message")

    # Print the progress bar
    cat(paste(c("\r  |",
                rep.int('=', current_step),
                rep.int(".", stepsLeft),
                sprintf("| Step %3d: %s", .global_step_counter, message),
                "\r"
                ), collapse=""), file = "")

    flush.console()
        # Reset the output redirection
    if (sink.number(type = "output") > 0) {
        sink(type = "output")
    }
    if (sink.number(type = "message") > 0) {
        sink(type = "message")
    }
}

.terminateProgressBar <- function() {
    .clearConsoleLine()
    .global_step_counter <<- 0
}

.initializeProgressBarOnLoad <- function(){
    # Global variable for progress bar step counter
    .global_step_counter <<- 0
}


###########################################
#  Databases versions printing functions  #
###########################################

.getCARDVersionAux <- function(){

    # Verify it has been downloaded and print version
    if(file.exists(paste(getOption("cardPath"),'/card-data/aro_index.tsv',sep=''))){

        # Read the JSON file as a text
        json_text <- readLines(paste(options("cardPath")[[1]],'/card-data/card.json',sep=''), warn = FALSE)

        # Combine the lines into a single string
        json_string <- paste(json_text, collapse = "")

        # Extract the _version value using regular expressions
        version <- sub('.*"_version":\\s*"([^"]*)".*', '\\1', json_string)

        # Extract the _timestamp value using regular expressions
        timestamp <- sub('.*"_timestamp":\\s*"([^"]*)".*', '\\1', json_string)
        date_part <- strsplit(timestamp, "T")[[1]][1]

        return(message(paste("CARD database version ", version, " (",date_part,")", sep="")))
    }
}

getCARDVersion <- function(){

    # Check if CARD is downloaded, do it if not
    if(!file.exists(paste(getOption("cardPath"),'/card-data/aro_index.tsv',sep=''))){
        updateCARDDataBase()
    } else {
        .getCARDVersionAux()
    }
}

getUniProtVersion <- function(){

    # URL to fetch the JSON data
    url <- "https://rest.uniprot.org/release-notes/search?query=*&sort=release_date+desc"

    # Use tryCatch to handle potential errors
    tryCatch({
        # Fetch the JSON content from the URL
        response <- httr::GET(
            url,
            httr::add_headers(
                `User-Agent` = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.3"
            )
        )

        # Check if the request was successful
        if (httr::status_code(response) == 200) {
        # Parse the JSON content
        json_data <- jsonlite::fromJSON(httr::content(response, "text", encoding = "UTF-8"), simplifyVector = FALSE)

        # Extract the first element in the results array
        first_result <- json_data$results[[2]]

        # Extract the title and releaseDate
        title <- first_result$title
        release_date <- first_result$releaseDate

        # Print the title and releaseDate
        message(paste(title, " (",release_date,")", sep=""))

        } else {
            message("Couldn't reach UniProt server, please try again.")
        }
    }, error = function(e) {
        message("Couldn't reach UniProt server, please try again.")
    })
}

getKEGGVersion <- function() {
    # URLs to fetch the HTML data
    url1 <- "https://www.genome.jp/kegg/docs/relnote.html"
    url2 <- "https://www.kegg.jp/kegg/rest/"

    # Function to fetch and parse the KEGG release version
    .fetchKEGGReleaseVersion <- function(url) {
        tryCatch({
            # Fetch the HTML content from the URL
            page <- httr::GET(url)

            # Check if the request was successful
            if (httr::status_code(page) == 200) {
                # Parse the HTML content using rvest
                content <- httr::content(page, as = "text", encoding = "UTF-8")
                html <- rvest::read_html(content)

                # Extract the release information
                release_info <- rvest::html_text(
                    rvest::html_node(
                        html,
                        xpath = "//h4[text()='Current release']/following-sibling::text()[1]"
                    ),
                    trim = TRUE
                )

                # Check if release information was found
                if (length(release_info) > 0) {
                    return(release_info)
                } else {
                    return(NULL)
                }
            } else {
                return(NULL)
            }
        }, error = function(e) {
            return(NULL)
        })
    }

    # Function to fetch and parse the KEGG API version
    .fetchKEGGAPIVersion <- function(url) {
        tryCatch({
            # Fetch the HTML content from the URL
            page <- httr::GET(url)

            # Check if the request was successful
            if (httr::status_code(page) == 200) {
                # Parse the HTML content using rvest
                content <- httr::content(page, as = "text", encoding = "UTF-8")
                html <- rvest::read_html(content)

                # Extract the update history information
                api_version_info <- rvest::html_text(
                    rvest::html_nodes(
                        html,
                        xpath = "//b[contains(text(), 'Update history')]/following-sibling::div/table/tr[1]/td[1]"
                    ),
                    trim = TRUE
                )

                # Check if API version information was found
                if (length(api_version_info) > 0) {
                    return(api_version_info)
                } else {
                    return(NULL)
                }
            } else {
                return(NULL)
            }
        }, error = function(e) {
            return(NULL)
        })
    }

    # Fetch and print the KEGG release version
    release_version <- .fetchKEGGReleaseVersion(url1)
    if (!is.null(release_version)) {
        message("KEGG database ", release_version)
    } else {
        message("Couldn't find the KEGG release information.")
    }

    # Fetch and print the KEGG API version
    api_version <- .fetchKEGGAPIVersion(url2)
    if (!is.null(api_version)) {
        message("KEGG API version ", api_version)
    } else {
        message("Couldn't find the KEGG API version information.")
    }
}


getNCBIVersion <- function() {
    # Helper function to extract and print the required fields
    .print_summary <- function(summary) {
        description <- summary[[2]]
        db_build <- summary[[4]]
        last_update <- sub(" .*", "", summary[[6]])  # Remove the hour part

        message('NCBI ', description, ' database ', db_build ,' (',last_update,')' , sep='')
    }

    # Fetch and print summaries for each database
    dbs <- c("protein", "gene", "nucleotide", "ipg")
    for (db in dbs) {
        summary <- entrez_db_summary(db)
        .print_summary(summary)
    }
}
