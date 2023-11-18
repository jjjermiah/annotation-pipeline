library(data.table)
library(AnnotationGx)
# inputfile <- "/home/bioinf/bhklab/jermiah/projects/annotationScripts/rawdata/ctrp/v20.meta.per_compound.txt"
# 
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

treatmentMetadata <- data.table::fread(treatmentMetadataFile, sep = "\t", check.names=TRUE)

# add prefix `CTRP.` to column names
# we will be using the CTRP.master_cpd_id as the primary key
names(treatmentMetadata) <- sapply(names(treatmentMetadata), function(x) paste0(toupper(dataset_name),".",x))

# clean the treatment names and create a column called `CTRP.cleanedTreatmentNames`
treatmentMetadata[, paste0("CTRP.cleanedTreatmentNames") := cleanCharacterStrings(CTRP.cpd_name)]

## ----------------- getPubChemCompound using all names ----------------- ##
message("running getPubChemCompound using all names")
compound_nameToCIDS <- 
    AnnotationGx::getPubChemCompound(
        treatmentMetadata[, CTRP.cleanedTreatmentNames], 
        from='name', 
        to='cids', 
        batch = FALSE,
        verbose = verbose,
        BPPARAM = BiocParallel::MulticoreParam(workers = THREADS, progressbar = TRUE, stop.on.error = FALSE)
    )


compound_successful_dt <- merge(
    treatmentMetadata, 
    compound_nameToCIDS, 
    by.x = "CTRP.cleanedTreatmentNames", 
    by.y = "name"
    )[,c("CTRP.master_cpd_id","cids")] # subset only the columns we want

failed_dt <- merge(
    treatmentMetadata, 
    data.table::rbindlist(attributes(compound_nameToCIDS)$failed, fill = TRUE),  # get all the failed queries
    by.x = "CTRP.cleanedTreatmentNames", 
    by.y = "query")

# ## ----------------- getPubChemSubstance using all failed names ----------------- ##
# message("running getPubChemSubstance using failed names")

# substancenameToCIDS <- 
#     AnnotationGx::getPubChemSubstance(
#         failed_dt$CTRP.cleanedTreatmentNames, 
#         from='name', 
#         to='sids', 
#         batch = FALSE)

# substancenameToCIDS <- substancenameToCIDS[!duplicated(name), ]
# message("running getPubChemSubstance using failed names")

# substance_successful_dt <- merge(
#     treatmentMetadata, 
#     substancenameToCIDS, 
#     by.x = "CTRP.cleanedTreatmentNames", 
#     by.y = "name"
#     )[,c("CTRP.master_cpd_id", "CID")]

# setnames(x=substance_successful_dt, old="CID", new="cids",skip_absent=TRUE)

# substance_failed_dt <- merge(
#     ctrp_treatment_meta, 
#     rbindlist(attributes(substancenameToCIDS)$failed, fill = TRUE), 
#     by.x = "CTRP.cleanedTreatmentNames", 
#     by.y = "query")
# message("running getPubChemSubstance using failed names")

## ----------------- getPubChemCompound using SMILES from failed ----------------- ##
message("running getPubChemCompound using SMILES from remaining failed.")
smilesToCIDS <- 
    AnnotationGx::getPubChemCompound(
        failed_dt$CTRP.cpd_smiles, 
        from='smiles', 
        to='cids', 
        batch = FALSE,
        verbose = verbose,
        BPPARAM = BiocParallel::MulticoreParam(workers = THREADS, progressbar = TRUE, stop.on.error = FALSE)
    )

smilesToCIDS <- smilesToCIDS[!duplicated(smiles), ]

smiles_successful_dt <- merge(
    treatmentMetadata, 
    smilesToCIDS, 
    by.x = "CTRP.cpd_smiles", 
    by.y = "smiles")

smiles_successful_dt <- smiles_successful_dt[cids != 0 & !duplicated(CTRP.master_cpd_id), c("CTRP.master_cpd_id", "cids")]


## ----------------- get list of failed ids ----------------- ##
failed_0 <- smiles_successful_dt[smiles_successful_dt$cids == 0, ] # some of these successful rows have cids = 0

if(length(attributes(smilesToCIDS)$failed) > 0){
    failed_1 <- 
        merge(
            treatmentMetadata, 
            data.table::rbindlist(attributes(smilesToCIDS)$failed, fill = TRUE), 
            by.x = "CTRP.cpd_smiles", 
            by.y = "query"
        )[,c("CTRP.master_cpd_id","CTRP.cpd_name")]
    failed <- data.table::rbindlist(list(failed_0, failed_1), fill = TRUE)[,c("CTRP.master_cpd_id","cids")]
}else failed <- failed_0[,c("CTRP.master_cpd_id","cids")]
    

## ----------------- compile lists----------------- ##
successful_CTRP.master_cpd_ids <- 
    data.table::rbindlist(
        list(
            compound_successful_dt, 
            smiles_successful_dt), 
        fill = TRUE
    )

failed <- failed[!duplicated(CTRP.master_cpd_id), ]

# create an object that is a list of the input data.table, successful, and failed
# this is to be saved as an RDS file
mappedCIDs <- list(treatmentMetadata, successful_CTRP.master_cpd_ids, failed)
names(mappedCIDs) <- c(
    paste0(dataset_name, ".treatment.metadata"), 
    paste0(dataset_name, ".successful.mapped.cids"), 
    paste0(dataset_name, ".failed.mapped.cids"))


# save the object as an RDS file
saveRDS(ctrp_cids, outputFile)

