## getPubchemAnnotations.R
library(data.table)
library(AnnotationGx)

# devtools::install_github("bhklab/AnnotationGx", ref = "jermiah", dependencies = TRUE, upgrade="never")
# get snakemake output
outputfile <- snakemake@output[['pubchem_annotation_db']]

# get the wildcard value
query <- snakemake@wildcards[['query']]

# get the params
# params <- snakemake@params[['query_type']]

# get snakemake number of threads
threads <- snakemake@threads

# BiocParallel::register(
#     BiocParallel::MulticoreParam(
#         workers=threads, 
#         progressbar = TRUE
#     )
# )
# AnnotationGx::getPubChemAnnotations
maxPages <- as.numeric(snakemake@params[['maxPages']])
if (is.na(maxPages) | maxPages == "") maxPages <- Inf
print(maxPages)

result <- 
    AnnotationGx::getPubChemAnnotations(
        header=query, 
        raw = as.logical(snakemake@params[['raw']]), 
        rawAnnotationDT = as.logical(snakemake@params[['rawAnnotationDT']]), 
        BPPARAM = BiocParallel::MulticoreParam(workers=threads, progressbar = TRUE),
        verbose = TRUE,
        maxPages = maxPages)

saveRDS(result, outputfile)
# headers <- AnnotationGx::getPubChemAnnotations()
# headers <- headers[Type == "Compound"]

# headers[Heading %like% 'ATC'          ``
#     | Heading %like% 'NSC' 
#     | Heading %like% 'CAS' 
#     | Heading %like% 'Drug Induced' 
#     | Heading %like% 'CTD' 
#     | Heading %like%  'FDA Approved' 
#     | Heading %like% 'Synonym'
#     | Heading %like% 'Clinical Trial']


# headers[Heading %like% 'ATC']
# CIDtoATC <- AnnotationGx::getPubChemAnnotations('ATC Code', raw = FALSE, rawAnnotationDT = TRUE)

# headers[Heading %like% 'NSC']
# CIDtoNSC <- AnnotationGx::getPubChemAnnotations('NSC Number', raw = FALSE, rawAnnotationDT = TRUE)

# headers[Heading %like% 'ChEMBL']
# CIDtoChEMBL <- AnnotationGx::getPubChemAnnotations('ChEMBL ID', raw = FALSE, rawAnnotationDT = TRUE)

# headers[Heading %like% 'Drug Induced']
# CIDtoDILI <- AnnotationGx::getPubChemAnnotations('Drug Induced Liver Injury',raw = FALSE, rawAnnotationDT = TRUE)

# headers[Heading %like% 'CAS']
# CIDtoCAS <- AnnotationGx::getPubChemAnnotations('CAS', raw = FALSE, rawAnnotationDT = TRUE)

# headers[Heading %like% 'FDA Approved Drugs']
# CIDtoFDA <- AnnotationGx::getPubChemAnnotations('FDA Approved Drugs', raw = FALSE, rawAnnotationDT = TRUE)

# headers[Heading %like% 'Syn']
# CIDtoNamesAndSyn <- AnnotationGx::getPubChemAnnotations('Names and Synonyms', raw = FALSE, rawAnnotationDT = TRUE)
# CIDtoSynAndIDs <- AnnotationGx::getPubChemAnnotations('Synonyms and Identifiers', raw = FALSE, rawAnnotationDT = TRUE)



# ## ---------- Parsing 
# ATC_anntns <- AnnotationGx:::.parseATCannotations(CIDtoATC)
# NSC_anntns <- AnnotationGx:::.parseNSCannotations(CIDtoNSC)

# DT <- CIDtoATC

# # sort DT on CID
# setkey(DT, CID)

