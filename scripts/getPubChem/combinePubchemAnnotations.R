## ------------------- Parse Snakemake Object ------------------- ##
if(exists("snakemake")){
    inputFiles <- snakemake@input[['pubchem_annotation_db']]
    inputPubchemProperties <- snakemake@input[['pubchem_properties']]
    annotationTypes <- sapply(inputFiles, function(x){
        file <- paste((basename(x)), collapse=" ")      # remove the directories from path
        file <- strsplit(file, "_")[[1]]                # split on underscore 
        strsplit(file[[length(file)]], "\\.")[[1]][[1]] # get last element and remove the extension
    })
    outputFile <- snakemake@output[['pubchemTreatmentAll']]
    save.image()
}

# read in each file and remove duplicated "cid" column, only take first two columns
data_ <- lapply(inputFiles, function(x){
    file <- readRDS(x)
    if(all(class(file)!='data.table')){             # sometimes the file is a list of data.tables
        suppressWarnings(
            file <- data.table::rbindlist(file, fill=TRUE)
        )    
    }
    file <- file[!duplicated(file$cid), 1:2]                      # remove duplicated "cid" column, only take first two columns
    return (file)
})

# merge all the data together on the "cid" column 
data <- Reduce(function(x, y) merge(x, y, by = "cid", all = TRUE), data_)
# rename "cid" to "CID"
data.table::setnames(data, "cid", "CID")

names(data) <- sapply(names(data), function(x) paste0("PubChem.", gsub(" ", "_",x) ))


inputPubchemProperties <- readRDS(inputPubchemProperties)

inputPubchemProperties$PubChem.Synonym <- lapply(inputPubchemProperties$PubChem.Synonym, as.list)



finalPubChemAnnotated <- data.table::as.data.table(merge(data, inputPubchemProperties, by = "PubChem.CID", all = TRUE))

data.table::fwrite(data, outputFile, sep = ",", quote = FALSE, na = "NA", row.names = FALSE, col.names = TRUE)

