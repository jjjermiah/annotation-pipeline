library(data.table)

library(httr)
library(jsonlite)
# possible list of bad characters:
#  - "(", ")", "[", "]", "{", "}", "<", ">", "/", "\\", "|", "?", "*", "+", ".", "^", "$"
# [\xb5], [], [ ,], [;], [:], [-], [+], [*], [%], [$], [#], [{], [}], [[], []], [|], [\\^], [/], [\\], [.], [_], [ ], [(], [)]
cleanCharacterStrings <- function(name){

    # # make sure name is a string
    name <- as.character(name)

    # if there is a colon like in "Cisplatin: 1 mg/mL (1.5 mM); 5 mM in DMSO"
    # remove everything after the colon
    name <- gsub(":.*", "", name)

    # remove ,  ;  -  +  *  $  %  #  ^  _  as well as any spaces
    name <- gsub("[\\,\\;\\-\\+\\*\\$\\%\\#\\^\\_\\s]", "", name, perl = TRUE)

    # remove hyphen 
    name <- gsub("-", "", name)

    # # remove substring of round brackets and contents
    name <- gsub("\\s*\\(.*\\)", "", name)

    # # remove substring of square brackets and contents
    name <- gsub("\\s*\\[.*\\]", "", name)

    # # remove substring of curly brackets and contents
    name <- gsub("\\s*\\{.*\\}", "", name)

    # # convert entire string to uppercase
    name <- toupper(name)

    print(name)

    # dealing with unicode characters 
    name <- gsub("Unicode", "", iconv(name, "LATIN1", "ASCII", "Unicode"), perl=TRUE)

    name
}



buildURL <- function(...) paste0(na.omit(unlist(list(...))), collapse='/')

parseJSON <- function(response, ..., encoding='UTF-8', query_only=FALSE) {
    if (isTRUE(query_only)) return(response)
    tryCatch({
        fromJSON(content(response, ..., as='text', type='JSON',
            encoding=encoding))
    },
    error=function(e) {
        fromJSON(content(response, ..., type='JSON', encoding=encoding))
    })
}

checkThrottlingStatus <- function(result, throttleMessage = FALSE){
    message <- headers(result)$`x-throttling-control`

    if (throttleMessage == TRUE){
        message(message)
    }
    matches <- regmatches(message, gregexpr("\\((.*?)%\\)", message))  # Extracts text within parentheses
    percentages <- gsub("\\(|%|\\)", "", unlist(matches[1:3]))
    # print(percentages)
    percentage <- max(as.numeric(percentages))
    if(as.integer(percentage) > 15 && as.integer(percentage) < 30){
        Sys.sleep(15)
    }else if (as.integer(percentage) > 30 && as.integer(percentage) < 50){
        Sys.sleep(20)
    }else if (as.integer(percentage) > 50 && as.integer(percentage) < 75) {
        print(paste0("Throttling at ", percentage, "%. Sleeping for 30 seconds."))
        Sys.sleep(30)
    }else if (as.integer(percentage) > 75) {
        print(paste0("Throttling at ", percentage, "%. Sleeping for 60 seconds."))
        Sys.sleep(60)
    }else{
        Sys.sleep(5)
    }   
}

getPubChemCHEMBL <- function(
    compound,
    url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound',
    output = 'JSON', 
    type = 'ChEMBL ID',
    timeout_s = 29,
    retries = 3,
    quiet = TRUE,
    throttleMessage = FALSE
    ){

        if(type == "DILI") queryURL <- paste0(buildURL(url, compound, output), '?heading=', "Drug Induced Liver Injury")
        else queryURL <- paste0(buildURL(url, compound, output), '?heading=', type)

        tryCatch({
            result <- RETRY('GET', URLencode(queryURL), times = retries, quiet = quiet)
        }, error=function(e) {
            print(paste0("Error: ", e$message))
            return(NULL)
        })

        checkThrottlingStatus(result, throttleMessage = throttleMessage)
        result <- parseJSON(result)

        if (type == 'ChEMBL ID') {
            result <- result$Record$Reference$SourceID
            result <- gsub("::Compound", "", result)
        }else if (type == 'NSC Number'){
            result <- result$Record$Reference$SourceID[1]
            result <- gsub(" ", "", result)
        }else if (type == 'DILI' || type =='Drug Induced Liver Injury'){
            if(length(result$Record$Section) == 0){
                result <- "NA"
                
            }else{
                dt_ <- as.data.table(result$Record$Section)
                dt_ <- as.data.table(dt_)$Section[[1]]
                dt_ <- as.data.table(dt_)$Section
                dt_ <- as.data.table(dt_)
                dt_ <- as.data.table(dt_)$Information
                # print(as.data.table(dt_)[1:3,  .(Name,unlist(Value))])
                section <- as.data.table(dt_)[1:3, "DILI" := paste0(unlist(Name), ":", unlist(Value))]

                # if any of the first 4 rows are NA, remove it 
                section <- section[!is.na(section$DILI)]


                section <- paste0(section[1:3, DILI], collapse= "; ")

                # create a list for each row as Name:Value string with no spaces and no new lines
                reference <- paste0("LTKBID:", result$Record$Reference$SourceID)
                result <- c(section, reference)
                result <- paste0(result, collapse = "; ")
            }
        }else if (type == 'CAS'){
            result <- result$Record$Reference$SourceID[1]
        }else if (type == 'ATC Code'){
            if(length(result$Record$Section) == 0){
                result <- "NA"
                
            }else{dt_ <- as.data.table(result$Record$Section)
            dt_ <- as.data.table(dt_)$Section[[1]]
            dt_ <- as.data.table(dt_)$Information
            dt_ <- as.data.table(dt_)$Value
            dt_ <- as.data.table(dt_[[1]])
            result <- paste0("ATC:", dt_$String)
            if(length(result) > 1){
                result <- paste0(result, collapse = "; ")
            }}
        }
        
        if (is.null(result)) result <- list(compound, "N/A")
        else result <- list(compound,result)
        names(result) <- c("cid", type)
        return(result)
    }


# # 'ATC Code'=return(.parseATCannotations(annotationDT)),
# # 'Drug Induced Liver Injury'=return(.parseDILIannotations(annotationDT)),
# # 'NSC Number'=return(.parseNSCannotations(annotationDT)),
# # 'CTD Chemical-Gene Interactions'=return(.parseCTDannotations(annotationDT)),
# # 'Names and Synonyms'=return(.parseNamesAndSynonyms(annotationDT)),
# # 'Synonyms and Identifiers'=return(.parseSynonymsAndIdentifiers(annotationDT)),
# # 'CAS'=return(.parseCASannotations(pageList)),


# # as.data.table(result$Record$Section[, c("TOCHeading", "Section")])[TOCHeading=="Names and Identifiers", Section][[1]][,'Section'][[3]][,'Information']


