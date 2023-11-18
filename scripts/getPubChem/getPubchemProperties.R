
## ------------------- Parse Snakemake Object ------------------- ##
if(exists("snakemake")){
    inputfile <- snakemake@input[['mappedCIDs']]
    outputFile <- snakemake@output[['pubchem_properties']]

    dataset_name <- snakemake@wildcards[['dataset']]           

    THREADS <- snakemake@threads
    snakemake@source("newfunctions.R")
    save.image()
}

data <- readRDS(inputfile)
mapped_cids <- data[[paste0(dataset_name, ".successful.mapped.cids")]][!is.na("cids")]

# TODO::ADD BPPARAM HERE
propertiesFromCID <- 
    AnnotationGx::getPubChemCompound(
        mapped_cids[, cids], 
        from='cid', 
        to='property', 
        properties=c('Title', 'MolecularFormula', 'InChIKey', 'CanonicalSMILES'))

CIDtoSynonyms <- 
    AnnotationGx::getPubChemCompound(
        mapped_cids[, cids], 
        from='cid', 
        to='synonyms')
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