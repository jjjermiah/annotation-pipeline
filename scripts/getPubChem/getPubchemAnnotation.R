library(data.table)

library(httr)
library(jsonlite)
## ------------------- Parse Snakemake Object ------------------- ## v
if(exists("snakemake")){
    inputfile <- snakemake@input[['mappedCIDs']]
    outputfile <- snakemake@output[['pubchem_annotation_dbs']]

    dataset_name <- snakemake@wildcards[['dataset']]           
    annotationType <- snakemake@wildcards[['annotationType']]
    THREADS <- snakemake@threads
    snakemake@source("newfunctions.R")
    save.image()
    
}

## -----------------  CHECK INPUTS ----------------- ##
# make sure inputfile and output file are strings 
stopifnot(is.character(inputfile))
stopifnot(is.character(outputfile))

annotations <- c('ChEMBL ID', 'NSC Number', 'Drug Induced Liver Injury', 'CAS', 'ATC Code')
stopifnot(is.character(annotationType))
stopifnot(annotationType %in% annotations)

## -----------------  READ INPUT DATA ----------------- ##
data <- readRDS(inputfile)

mapped_cids <- data[[paste0(dataset_name, ".successful.mapped.cids")]][!is.na("cids")]

## -----------------  GET ANNOTATIONS  ----------------- ##
## should return a list of data.tables
## including the cid to merge back on 
## and the annotation
result <- 
    BiocParallel::bptry(
        BiocParallel::bplapply(
                mapped_cids[, cids], 
                function(CID) {as.data.table(AnnotationGx::getPubChemAnnotation(CID, header = annotationType))},
                BPPARAM = BiocParallel::MulticoreParam(workers = THREADS, progressbar = TRUE, stop.on.error = FALSE)
        )
    )

# Get the successful runs
successful <- result[BiocParallel::bpok(result)]
failed <- result[!BiocParallel::bpok(result)] # unused


CIDS_to_Annotation <- rbindlist(successful, fill = TRUE)

# # save result to RDS file 
print(paste0("Saving result to ", outputfile))
saveRDS(CIDS_to_Annotation, outputfile)
