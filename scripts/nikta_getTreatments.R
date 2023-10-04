library(AnnotationGx)
library(BiocParallel)
library(data.table)
library(httr)
library(stringr)

rbind_annot <- readRDS("/home/bioinf/bhklab/jermiah/projects/annotationScripts/rawdata/rbind_annot.rds") 
rbind_annot <- as.data.table(rbind_annot)


# 4. Finding ChembleIDs from cids
processQuery <- function(compound) {
  
  url <- 'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound'
  output <- 'JSON'
  type <- 'ChEMBL+ID'
  queryURL <- paste0(.buildURL(url, compound, output),
                     '?heading=', type)
  encodedQueryURL <- URLencode(queryURL)
  queryRes <- RETRY('GET', encodedQueryURL, timeout(29), times = 10, quiet = FALSE)
  AnnotationGx:::.checkThrottlingStatus(queryRes)
  result <- parseJSON(queryRes, as = 'text')
  reference <- result$Record$Reference
  sourceID <- reference$SourceID


  # if sourceID is empty, then set it to NA
  if (length(sourceID) == 0) desiredString <- NA
  else desiredString <- gsub("::Compound", "", sourceID[1])

  result_df <- list(cids = compound, ChEMBL_ID= desiredString, query = encodedQueryURL)
  
  return(as.data.table(result_df,fill=TRUE))
  
}

# To avoid the "Error in the HTTP2 framing layer"
httr::set_config(httr::config(http_version <- 0))
chemble_ids <- rbindlist(lapply(rbind_annot$ref_cid, processQuery))

# subset the rows where ChEMBL_ID is not NA
successful_ChEMBL <- chemble_ids[!is.na(ChEMBL_ID), c("cids", "ChEMBL_ID")]
failed_ChEMBL <- chemble_ids[is.na(ChEMBL_ID)] 
failed_ChEMBL

# Adding Chemble ids to the drugnames
chemble_cid = merge( x= rbind_annot, y = chemble_ids, by.x = "ref_cid", by.y = "cids",)

# Got the below line working using "https://github.com/ropensci/rtweet/issues/229"
# According to which I had to change my DNS to 8.8.8.8
all_moas = getChemblAllMechanisms()

final_annot = merge(all_moas, chemble_cid, by.x = "molecule_chembl_id", by.y = "ChEMBL_ID") %>%
              subset(select = c("org_drug", 
                                "ref_title",
                                "ref_cid", 
                                "molecule_chembl_id", 
                                "parent_molecule_chembl_id",
                                "action_type", 
                                "mechanism_of_action", 
                                "molecular_mechanism",
                                "mechanism_comment"))

# Save the results
saveRDS(final_annot, "~/Documents/PHD/Drug response project/drug_screen/Drug annotation/final_annot.rds")