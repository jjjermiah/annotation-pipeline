library(data.table)
library(AnnotationGx)
# inputfile <- "/home/bioinf/bhklab/jermiah/projects/annotationScripts/rawdata/ctrp/v20.meta.per_compound.txt"
# 

# get snakemake input and output files
inputfile <- snakemake@input[['dataset_annotation_file']]
outputfile <- snakemake@output[['treatment_Pubchem_data']]
dataset_name <- snakemake@wildcards[['dataset']]

ctrp_treatment_meta <- as.data.table(read.delim(inputfile, header=TRUE, sep="\t"))
# Check for matches differing only by case, spaces, or the following characters:
# possible list of bad characters:
#  - "(", ")", "[", "]", "{", "}", "<", ">", "/", "\\", "|", "?", "*", "+", ".", "^", "$"
# [\xb5], [], [ ,], [;], [:], [-], [+], [*], [%], [$], [#], [{], [}], [[], []], [|], [\\^], [/], [\\], [.], [_], [ ], [(], [)]

# remove all characters in the list above
cleanTreatmentNames <- function(name){
    # make sure name is a string
    name <- as.character(name)
    # convert entire string to uppercase
    name <- toupper(name)

    # remove substring of round brackets and contents
    name <- gsub("\\s*\\(.*\\)", "", name)

    # remove substring of square brackets and contents
    name <- gsub("\\s*\\[.*\\]", "", name)

    # remove substring of curly brackets and contents
    name <- gsub("\\s*\\{.*\\}", "", name)

    # # remove colon and replace with no space
    # name <- gsub(":", "", name)

    # # remove hyphen and replace with no space
    # name <- gsub("-", "", name)

    # # remove plus and replace with no space
    # name <- gsub("\\+", "", name)

    # # remove asterisk and replace with no space
    # name <- gsub("\\*", "", name)

    # remove all spaces
    # name <- gsub("\\s", "", name)

    name
}

treatments <- as.data.table(sapply(ctrp_treatment_meta$cpd_name, cleanTreatmentNames),keep.rownames=T)
names(treatments) <- c("cpd_name", "treatment_name_cleaned")
ctrp_treatment_meta <- merge(ctrp_treatment_meta, treatments, by = "cpd_name")

# BiocParallel::register(BiocParallel::SerialParam())
BiocParallel::register(BiocParallel::MulticoreParam(workers=1))


## ----------------- getPubChemCompound using all names ----------------- ##
print("running getPubChemCompound using all names")
compound_nameToCIDS <- 
    AnnotationGx::getPubChemCompound(
        ctrp_treatment_meta$treatment_name_cleaned, 
        from='name', 
        to='cids', 
        batch = FALSE)

compound_nameToCIDS <- compound_nameToCIDS[!duplicated(name), ]

compound_successful_dt <- merge(
    ctrp_treatment_meta, 
    compound_nameToCIDS, 
    by.x = "treatment_name_cleaned", 
    by.y = "name"
    )[,c("master_cpd_id","cids")] # subset only the columns we want

failed_dt <- merge(
    ctrp_treatment_meta, 
    rbindlist(attributes(compound_nameToCIDS)$failed, fill = TRUE),  # get all the failed queries
    by.x = "treatment_name_cleaned", 
    by.y = "query")

## ----------------- getPubChemSubstance using all failed names ----------------- ##
substancenameToCIDS <- 
    AnnotationGx::getPubChemSubstance(
        failed_dt$treatment_name_cleaned, 
        from='name', 
        to='cids', 
        batch = FALSE)

substancenameToCIDS <- substancenameToCIDS[!duplicated(name), ]

substance_successful_dt <- merge(
    ctrp_treatment_meta, 
    substancenameToCIDS, 
    by.x = "treatment_name_cleaned", 
    by.y = "name"
    )[,c("master_cpd_id", "CID")]

setnames(x=substance_successful_dt, old="CID", new="cids",skip_absent=TRUE)

substance_failed_dt <- merge(
    ctrp_treatment_meta, 
    rbindlist(attributes(substancenameToCIDS)$failed, fill = TRUE), 
    by.x = "treatment_name_cleaned", 
    by.y = "query")

## ----------------- getPubChemCompound using SMILES from failed ----------------- ##
smilesToCIDS <- 
    AnnotationGx::getPubChemCompound(
        substance_failed_dt$cpd_smiles, 
        from='smiles', 
        to='cids', 
        batch = FALSE)

smilesToCIDS <- smilesToCIDS[!duplicated(smiles), ]

smiles_successful_dt <- merge(
    ctrp_treatment_meta, 
    smilesToCIDS, 
    by.x = "cpd_smiles", 
    by.y = "smiles")

smiles_successful_dt <- smiles_successful_dt[cids != 0 & !duplicated(master_cpd_id), c("master_cpd_id", "cids")]

failed_0 <- smiles_successful_dt[smiles_successful_dt$cids == 0, ] # some of these successful rows have cids = 0
failed_1 <- merge(ctrp_treatment_meta, rbindlist(attributes(smilesToCIDS)$failed, fill = TRUE), by.x = "cpd_smiles", by.y = "query")[,c("master_cpd_id","cpd_name")]

## ----------------- get list of failed ids ----------------- ##
failed_master_cpd_ids <- rbindlist(list(failed_0, failed_1), fill = TRUE)[,c("master_cpd_id","cpd_name")]

## ----------------- compile lists----------------- ##
successful_master_cpd_ids <- 
    rbindlist(
        list(
            compound_successful_dt, 
            substance_successful_dt, 
            smiles_successful_dt), 
        fill = TRUE
    )[!duplicated(master_cpd_id),c("master_cpd_id","cids")]

failed_master_cpd_ids <- failed_master_cpd_ids[!duplicated(master_cpd_id), ]

# create an object that is a list of the input data.table, successful, and failed
# this is to be saved as an RDS file
ctrp_cids <- list(
    ctrp.treatment.metadata = ctrp_treatment_meta,
    successful_master_cpd_ids = successful_master_cpd_ids,
    failed_master_cpd_ids = failed_master_cpd_ids
)

# save the object as an RDS file
saveRDS(ctrp_cids, outputfile)