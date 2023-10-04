library(data.table)
library(AnnotationGx)
library(httr)
inputfile <- "/home/bioinf/bhklab/jermiah/projects/annotationScripts/results/ctrp/pubchem_preproccessed.RDS"

data <- readRDS(inputfile)

successful_cids <- data$successful_master_cpd_ids
names(successful_cids) <- c("ctrp.treatmentID", "cids")
successful_cids <- successful_cids[cids != "NULL"]
successful_cids <- successful_cids[, cids := unlist(cids)]

BiocParallel::register(BiocParallel::MulticoreParam(workers=8, progressbar = TRUE))
# get the other properties for the successful cids
propertiesFromCID <- 
    AnnotationGx::getPubChemCompound(
        successful_cids$cids[1:10], 
        from='cid', 
        to='property', 
        properties=c('Title', 'MolecularFormula', 'InChI', 'InChIKey', 'CanonicalSMILES')
    )
# reorder columns to be 'CID', 'Title', 'MolecularFormula', 'InChI', 'InChIKey', 'CanonicalSMILES'
propertiesFromCID <- propertiesFromCID[,c('CID', 'Title', 'MolecularFormula', 'InChI', 'InChIKey', 'CanonicalSMILES')]


getPubChemPugView <- 
    function(
        compound,
        header = 'data', 
        type = 'ChEMBL+ID', 
        output = 'JSON',
        raw = FALSE
        ) {
    if (header == 'Available') {
        queryURL <-
            'https://pubchem.ncbi.nlm.nih.gov/rest/pug/annotations/headings/JSON'
    } else if (header == 'data') {
        url <- 'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound'
        queryURL <- paste0(.buildURL(url, compound, output),
            '?heading_type=', type)
    }

    encodedQueryURL <- URLencode(queryURL)

    queryRes <- httr::RETRY('GET', encodedQueryURL, timeout(29), times = 10, quiet = TRUE)
    AnnotationGx:::.checkThrottlingStatus(queryRes)

    if (raw) return(queryRes)

    
    tryCatch({
        result <- parseJSON(queryRes, as = 'text')
        reference <- as.data.table(result$Record$Reference)

        # subset reference data.table where "SourceName" is "ChEMBL" and "SourceID" contains "CHEMBL"
        ChEMBL_ID <- reference[SourceName=="ChEMBL" & grepl("CHEMBL", SourceID), SourceID]
        # if sourceID is empty, then set it to NA
        if (length(ChEMBL_ID) == 0) desiredString <- NA
        else desiredString <- gsub("::Compound", "", ChEMBL_ID[1])
        result_df <- list(cids = compound, ChEMBL_ID = desiredString, query = encodedQueryURL)
        return(as.data.table(result_df,fill=TRUE))

    }, error = function(e) {
        result <- NA
        result_df <- list(cids = compound, ChEMBL_ID = NA, query = as.character(e))
        return(as.data.table(result_df,fill=TRUE))
    })
}

# getPubChemPugView(compound = 132971)

queryList <- 
    BiocParallel::bplapply(successful_cids[cids != "NULL", cids], function(x){
        getPubChemPugView(
            compound = as.character(x), 
            header = 'data', 
            type = 'ChEMBL+ID', 
            output = 'JSON', 
            raw = FALSE)
        }
    )

queryRes <- rbindlist(queryList, fill = TRUE)
# convert the "cids" column to character strings
queryRes <- queryRes[, cids := as.integer(cids)]

merge(queryRes, successful_cids, by = "cids")




# 'ATC Code'=return(.parseATCannotations(annotationDT)),
# 'Drug Induced Liver Injury'=return(.parseDILIannotations(annotationDT)),
# 'NSC Number'=return(.parseNSCannotations(annotationDT)),
# 'CTD Chemical-Gene Interactions'=return(.parseCTDannotations(annotationDT)),
# 'Names and Synonyms'=return(.parseNamesAndSynonyms(annotationDT)),
# 'Synonyms and Identifiers'=return(.parseSynonymsAndIdentifiers(annotationDT)),
# 'CAS'=return(.parseCASannotations(pageList)),


# as.data.table(result$Record$Section[, c("TOCHeading", "Section")])[TOCHeading=="Names and Identifiers", Section][[1]][,'Section'][[3]][,'Information']


