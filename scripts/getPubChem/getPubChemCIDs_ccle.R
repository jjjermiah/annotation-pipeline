## ------------------- Parse Snakemake Object ------------------- ##
if(exists("snakemake")){
    treatmentMetadataFile <- snakemake@input[['dataset_treatment_metadata_file']]
    dataset_name <- snakemake@wildcards[['dataset']]            # dataset name
    outputFile <- snakemake@output[['mappedCIDs']]

    THREADS <- snakemake@threads
    verbose <- F

    snakemake@source("newfunctions.R")

    save.image()
}

treatmentMetadata <- data.table::fread(treatmentMetadataFile, header = TRUE,  encoding = "Latin-1") # need to use Latin-1 because of weird characters after Panobinostat

treatmentMetadata <-
    treatmentMetadata[,
        .(  CCLE.treatmentID = `Compound (code or generic name)`,
            CCLE.target = `Target(s)`,
            CCLE.mechanismOfAction = `Mechanism of action`,
            CCLE.class = `Class`,
            CCLE.highestClinicalTrialPhase = `Highest Phase`,
            CCLE.treatmentSourceOrganization = `Organization`)]

treatmentMetadata[, CCLE.cleanedTreatmentName := cleanCharacterStrings(CCLE.treatmentID)]

message("running getPubChemCompound using all names")
compound_nameToCIDS <-
    AnnotationGx::getPubChemCompound(
        treatmentMetadata[, CCLE.cleanedTreatmentName],
        from='name',
        to='cids',
        batch = FALSE,
        verbose = verbose,
        BPPARAM = BiocParallel::MulticoreParam(workers = THREADS, progressbar = TRUE, stop.on.error = FALSE)
    )
compound_nameToCIDS <- compound_nameToCIDS[!duplicated(name), ]

failed <- attributes(compound_nameToCIDS)$failed # note: as of writing, CCLE has no failed queries. 

# create a list of the input data.table, successful, and failed
mappedCIDs <- list(treatmentMetadata, compound_nameToCIDS, failed)
names(mappedCIDs) <- c(
    paste0(dataset_name, ".treatment.metadata"), 
    paste0(dataset_name, ".successful.mapped.cids"), 
    paste0(dataset_name, ".failed.mapped.cids"))

saveRDS(mappedCIDs, outputFile)