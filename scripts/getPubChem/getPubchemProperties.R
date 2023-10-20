
## ------------------- Parse Snakemake Object ------------------- ##
if(exists("snakemake")){
    input <- snakemake@input[['treatment_Pubchem_data']]

    outputFile <- snakemake@output[['pubchem_properties']]
}else{
    input <- "/home/bioinf/bhklab/jermiah/projects/annotationScripts/procdata/ctrp/pubchem_preproccessed.RDS"
}

data <- readRDS(input)

CIDS <- 
    unlist(
        lapply(
            data$successful_master_cpd_ids$cids, 
            function(x) ifelse(is.null(x) || is.null(unlist(x)), NA, x)))
CIDS <- CIDS

# TODO::ADD BPPARAM HERE
propertiesFromCID <- 
    AnnotationGx::getPubChemCompound(
        CIDS, 
        from='cid', 
        to='property', 
        properties=c('Title', 'MolecularFormula', 'InChIKey', 'CanonicalSMILES'))


CIDtoSynonyms <- AnnotationGx::getPubChemCompound(CIDS, from='cid', to='synonyms')
# CIDtoSIDs <- AnnotationGx::getPubChemCompound(CIDS, from='cid', to='sids')

# merge synonyms and properties
propertiesFromCID <- merge(propertiesFromCID, CIDtoSynonyms, by.x='CID', by.y='CID')

# add prefix `PubChem.` to column names
names(propertiesFromCID) <- sapply(names(propertiesFromCID), function(x) paste0("PubChem.",x))

# set column order 
propertiesFromCID <- 
    propertiesFromCID[,c('PubChem.CID', 'PubChem.Title', 'PubChem.MolecularFormula', 'PubChem.InChIKey', 'PubChem.CanonicalSMILES', 'PubChem.Synonym')]

saveRDS(propertiesFromCID, outputFile)

# save.image()