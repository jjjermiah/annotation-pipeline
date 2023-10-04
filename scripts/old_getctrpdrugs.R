library(data.table)
library(tidyverse)
serialRun <- BiocParallel::SerialParam()
parallelRun<- BiocParallel::MulticoreParam(workers=4)
BiocParallel::register(serialRun)
## PUBCHEM

# going to test using the drugs from CTRPv2
# RAWDATA
ctrp_treatment_meta <- as.data.table(read.delim("rawdata/CTRPv2/v20.meta.per_compound.txt", header=TRUE, sep="\t"))
x <- ctrp_treatment_meta$cpd_name
nameToCIDS <- AnnotationGx::getPubChemCompound(x, from='name', to='cids')
nameToCIDS <- nameToCIDS[!duplicated(name), ]
nameToCIDS


# Get the failed queries, unnest the failure information to view potential reasons for failure
dt <- rbindlist(attributes(nameToCIDS)$failed, fill = TRUE)
failedq_dt <- as.data.table(unnest_wider(dt, failure))
failedq_dt


########################################################
####### Retrying for failed queries
# get the elements of dt where there is a colon
dt_colon <- failedq_dt[grepl(":", failedq_dt$query), query]

# format similar to "JQ-1:UNC0638 (2:1 mol/mol)", "carboplatin:UNC0638 (2:1 mol/mol)", 
#       "vorinostat:carboplatin (1:1 mol/mol)", "serdemetan:SCH-529074 (1:1 mol/mol)"   
# returns a list of elements after parsing each row dt_colon
processed_list <- unique(unlist(
    lapply(dt_colon, function(row) {
        # Step 1: Remove (#:# mol/mol) pattern
        row <- gsub("\\s*\\(\\d+:\\d+ mol/mol\\)", "", row)

        # Step 2: Separate string on ":" and append to a list
        elements <- strsplit(row, ":")[[1]]

        return(elements)
        }
    )
))
# Many of the combinations are two drugs, where the individual drugs are successful queries in nameToCIDS
# in get the elements of processed_list that are NOT in nameToCIDS$name
missing_drugs <- processed_list[!processed_list %in% nameToCIDS$name]

# elements from dt that dont have colons
dt_No_colon <- dt[!grepl(":", dt$query), query]

# combine both lists to try again using the smiles (provided by CTRP metadata) this time 
missing_drugs <- c(missing_drugs, dt_No_colon)
missing_drugs_smiles <- ctrp_treatment_meta[cpd_name %in% missing_drugs, cpd_smiles]
smilesToCIDS <- AnnotationGx::getPubChemCompound(missing_drugs_smiles, from='smiles', to='cids')
smilesToCIDS <- smilesToCIDS[!duplicated(smiles), ]

# get the failed queries, unnest the failure information to view potential reasons for failure
failedqueries<- attributes(smilesToCIDS)$failed
dt <- rbindlist(failedqueries, fill = TRUE)
failedq_smiles_dt <- as.data.table(unnest_wider(dt, failure))
failedq_smiles_dt

# get the rest of the metadata for the remaining failed drugs
remaining_failed_dt <- merge(failedq_smiles_dt, ctrp_treatment_meta, by.x='query', by.y='cpd_smiles')
remaining_failed_dt

########################################################
#### Working with the successful CIDS

##### data.table of the retrieved CIDS
successfulSmilesToCIDS <-  merge(smilesToCIDS, ctrp_treatment_meta, by.x='smiles', by.y='cpd_smiles')
successfulSmilesToCIDS <- successfulSmilesToCIDS[,c('cpd_name', 'cids')]
setnames(successfulSmilesToCIDS, 'cpd_name', 'name')
successfulCIDs <- rbind(nameToCIDS,successfulSmilesToCIDS)

# Use the successfulCIDs to get the other properties
propertiesFromCID <- AnnotationGx::getPubChemCompound(successfulCIDs$cids, from='cid', to='property', properties=c('Title', 'MolecularFormula', 'InChI', 'InChIKey', 'CanonicalSMILES'))
failedqueries <-  attributes(propertiesFromCID)$failed$`2`$query
failed_retry_propertiesFromCID <- AnnotationGx::getPubChemCompound(failedqueries, from='cid', to='property', properties=c('Title', 'MolecularFormula', 'InChI', 'InChIKey', 'CanonicalSMILES'), batch=FALSE)

# drop the 'cid' column from failed_propertiesFromCID
# rbind the successful properties with the retried failed properties
successfulproperties <- rbind(propertiesFromCID, failed_retry_propertiesFromCID[,c('CID', 'Title', 'MolecularFormula', 'InChI', 'InChIKey', 'CanonicalSMILES')])

# also get the SIDS (Substance Identifier from PubChem) for each CID
CIDtoSIDs <- AnnotationGx::getPubChemCompound(successfulproperties$CID, from='cid', to='sids')
CIDtoSynonyms <- AnnotationGx::getPubChemCompound(successfulproperties$CID, from='cid', to='synonyms')

# merge the properties with the successfulCIDs on CID
successfulproperties <- merge(successfulproperties, CIDtoSIDs)
propertiesFromCID <- merge(successfulCIDs, successfulproperties, by.x='cids', by.y='CID')

#reorder propertiesFromCID to have CIDS, Title, MolecularFormula, InChI, InChIKey, CanonicalSMILES
propertiesFromCID <- propertiesFromCID[,c('name','Title', 'cids',  'SID', 'MolecularFormula',  'InChIKey', 'CanonicalSMILES')]

propertiesFromCID


#########################################################
#### Get other annotations 
# https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/65064/JSON?heading=ChEMBL+ID
headers <- AnnotationGx::getPubChemAnnotations()
headers[Heading %like% 'ATC' 
    | Heading %like% 'NSC' 
    | Heading %like% 'CAS' 
    | Heading %like% 'Drug Induced' 
    | Heading %like% 'CTD' 
    | Heading %like%  'FDA Approved' 
    | Heading %like% 'Synonym'
    | Heading %like% 'Clinical Trial']
#
successfulCIDs$cids[1:10]
compound <- successfulCIDs$cids[1]
url <- 'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound'
output <- 'JSON'
type <- 'ChEMBL+ID'

library(httr)
library(AnnotationGx)

queryURL <- paste0(AnnotationGx::.buildURL(url, compound, output))
queryURL <- paste0(AnnotationGx::.buildURL(url, compound, output), '?heading=', type)
encodedQueryURL <- URLencode(queryURL)
result <- RETRY('GET', encodedQueryURL, timeout(29), times = 3, quiet = TRUE)
AnnotationGx:::.checkThrottlingStatus(result)
AnnotationGx::parseJSON(result)[[1]][[1]]

# Define the function to process each compound
processQuery <- function(compound, type) {
  queryURL <- paste0(.buildURL(url, compound, output),
                     '?heading=', type)
  encodedQueryURL <- URLencode(queryURL)
  queryRes <- RETRY('GET', encodedQueryURL, timeout(29), times = 3, quiet = TRUE)
  result <- parseJSON(queryRes, as = 'text')

  if (type=="ChEMBL+ID"){
    reference <- result$Record$Reference
    sourceID <- reference$SourceID
    desiredString <- gsub("::Compound", "", sourceID[1])
    result_df <- data.table(cids = compound, ChEMBL_ID= desiredString)
  }


  return(result_df)
}

parallelRun<- BiocParallel::MulticoreParam(workers=4)
BiocParallel::register(parallelRun)
# Getting CHEMBL IDs, may take some time.
CIDtoCHEMBL <- rbindlist(BiocParallel::bplapply(successfulCIDs$cids, processQuery, type = 'ChEMBL+ID'))
missing_ChEMBLIDS <- sum(is.na(CIDtoCHEMBL$ChEMBL_ID))
# final_dt <- merge(CIDtoCHEMBL, propertiesFromCID, by='cids')

# GET THE CHEMBL MECHANISMS(result <- getChemblAllMechanisms())
(result <- getChemblAllMechanisms())
result <- unique(result[,c('ChEMBL_ID',  "parent_molecule_chembl_id","action_type", "mechanism_of_action", "molecular_mechanism","mechanism_comment")])
result <-  result[!duplicated(ChEMBL_ID)]
CIDtoCHEMBL <-  merge(CIDtoCHEMBL, result, by = "ChEMBL_ID", all.x=TRUE)
cols <- c("cids", 
          "ChEMBL_ID", 
          "parent_molecule_chembl_id",
          "action_type", 
          "mechanism_of_action", 
          "molecular_mechanism",
          "mechanism_comment")
CIDtoCHEMBL[, ..cols]


final_ctrp_drug_dt <- merge(propertiesFromCID, CIDtoCHEMBL,by = "cids", all=TRUE)


# Convert all columns to character strings
final_ctrp_drug_dt[, names(final_ctrp_drug_dt) := lapply(.SD, as.character), .SDcols = names(final_ctrp_drug_dt)]

# Retry fwrite()
data.table::fwrite(final_ctrp_drug_dt, 'ctrp_drug_metadata.csv')

# CIDtoCHEMBL1 <- rbindlist(BiocParallel::bplapply(successfulCIDs$cids[1:100], processQuery, type = 'ChEMBL+ID'))
# CIDtoCHEMBL2 <- rbindlist(BiocParallel::bplapply(successfulCIDs$cids[101:200], processQuery, type = 'ChEMBL+ID'))
# CIDtoCHEMBL3 <- rbindlist(BiocParallel::bplapply(successfulCIDs$cids[201:300], processQuery, type = 'ChEMBL+ID'))
# CIDtoCHEMBL4 <- rbindlist(BiocParallel::bplapply(successfulCIDs$cids[301:400], processQuery, type = 'ChEMBL+ID'))
# CIDtoCHEMBL5 <- rbindlist(BiocParallel::bplapply(successfulCIDs$cids[401:length(successfulCIDs$cids)], processQuery, type = 'ChEMBL+ID'))

cols <- c("name", 
          "cids", 
          "ChEMBL_ID", 
          "parent_molecule_chembl_id",
          "action_type", 
          "mechanism_of_action", 
          "molecular_mechanism",
          "mechanism_comment")

# data.table::fwrite(CIDtoCHEMBL1, 'ctrp_cid_chemb
########################################################
#### rcellminer useful tools
pak::pkg_install("CBIIT/rcellminerUtilsCDB")
dbName<-'ctrp'
availableInDb <- !is.na(rcellminerUtilsCDB::drugSynonymTab[, dbName])
availableInDb
 dbDrugSynTab <- rcellminerUtilsCDB::drugSynonymTab[availableInDb, ]
dbDrugSynTab
dbDrugSynTab[,ctrp]
dbDrugSynTab[,'ctrp']
as.data.table(dbDrugSynTab[,'ctrp'])
dbDrugSynTab