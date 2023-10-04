#  GetCellosaurus_ctrp.R 
library(tidyverse)
library(data.table)
library(AnnotationGx)
library(Cellosaurus)


########
# LOAD FILES
ctrpcells <- "/home/bioinf/bhklab/jermiah/projects/annotationScripts/rawdata/CTRPv2/v20.meta.per_cell_line.txt"
ctrpcells <- fread(ctrpcells, header = T, sep = "\t")

# Use Cellosaurus functions from new package:  https://r.acidgenomics.com/packages/cellosaurus/articles/Cellosaurus.html
object <- Cellosaurus::Cellosaurus()

# get all cellosaurusdata into data.table
cellosaurusData <- as.data.table(object)

# use mapping function, resulting ctrp_mapped is a named 'character' converted to a data.table

ctrp_mapped <- as.data.table(Cellosaurus::mapCells(object, cells = ctrpcells$ccl_name),keep.rownames = T)
names(ctrp_mapped) <- c("ctrp.sampleid", "cellosaurus.accession.ID") # set the names of the columns


########
# merge on cellosaurusData$accession and ctrp_mapped$rn
# cellosaurusData$accession is the accession number of the cell line
# ctrp_mapped$cellosaurus.accession.ID is the rowname of the ctrp_mapped data.table
merged_dt <- merge(cellosaurusData, ctrp_mapped, by.x = "accession", by.y = "cellosaurus.accession.ID", all.y = TRUE)

## pick columns that we want to keep
# names(merged_dt)
#  [1] "ctrp.sampleid" "cellLineName"   "synonyms"       "accession"     
#  [5] "diseases"       "samplingSite"   "oncotreeName"   "oncotreeTissue"
#  [9] "isCancer"       "category"       "depmapId"       "atccId"        
# [13] "sexOfCell"      "ageAtSampling" 

# list of all the columns to exclude from resulting data.table
# [1] "comments"                    "crossReferences"            
#  [3] "date"                        "hierarchy"                  
#  [5] "isContaminated"              "isProblematic"              
#  [7] "misspellings"                "msiStatus"                  
#  [9] "ncbiTaxonomyId"              "ncitDiseaseId"              
# [11] "ncitDiseaseName"             "oncotreeCode"               
# [13] "oncotreeLevel"               "oncotreeMainType"           
# [15] "oncotreeParent"              "organism"                   
# [17] "originateFromSameIndividual" "population"                 
# [19] "referencesIdentifiers"       "sangerModelId"              
# [21] "secondaryAccession"          "strProfileData"             
# [23] "webPages"  

results_metadata_dt <- merged_dt[, c(
    "ctrp.sampleid",  "cellLineName", "synonyms", "accession", 
    "diseases", "samplingSite", "oncotreeName", "oncotreeTissue", "isCancer", 
    "category", "depmapId", "atccId", "sexOfCell", "ageAtSampling"
    )]

results <- as.data.table(tidyr::unnest_wider(results_metadata_dt, diseases))


##################################################################
########### Using AnnotationGx


#### USING THE CTRP PROVIDED CELL_LINE NAMES
agx_results <- 
    AnnotationGx::getCellosaurusAPI(
        cl_names = standardizeCells(ctrpcells$ccl_name), 
        fields = c("ac", "sy", "misspelling", "din", "ca", "sx", "ag", "sampling-site", "metastatic-site"),
        GETfxn = "search/cell-line?",
        querydomain = "id:"
    )
cleaned <- AnnotationGx::cleanCellosaurusResponse(agx_results,GETfxn = "search/cell-line?")

##### USING THE MAPPED CELLOSUARUS ACCESSIONS FROM THE CTRP PROVIDED CELL_LINE NAMES
agx_results <- 
    AnnotationGx::getCellosaurusAPI(
        cl_names = ctrp_mapped$cellosaurus.accession.ID[1:10], 
        fields = c("ac", "sy", "misspelling", "din", "ca", "sx", "ag", "sampling-site", "metastatic-site"),
        GETfxn = "search/cell-line?",
        querydomain = "ac:"
    )
cleaned <- AnnotationGx::cleanCellosaurusResponse(agx_results, GETfxn="search/cell-line?")

##### USING THE MAPPED CELLOSUARUS ACCESSIONS FROM THE CTRP PROVIDED CELL_LINE NAMES AND DIFFERENT GETFXN
# cell-line/
responseList <- getCellosaurusAPI(ctrp_mapped$cellosaurus.accession.ID[1:10], GETfxn="cell-line/")
cleaned <- cleanCellosaurusResponse(responseList, GETfxn="cell-line/")
cleaned

# TEST TO SEE IF ANNOTATIONGX ANNOTATIONS AND CELLOSAURUS PACKAGE ANNOTATIONS ARE SIMILAR 

# sort both data.tables by cellLineName
cleaned <- cleaned[order(cellLine)]
results <- results[order(ctrp.sampleid)]


# subset the results data.table using only the ctrp.sampleid that are in the cleaned data.table cellLine column
results_subset <- results[ctrp.sampleid %in% cleaned$cellLine]

# join results_subset and cleaned on cellLine and only keep results_subset$accession and cleaned$AC
joined <- merge(results_subset, cleaned, by.x = "ctrp.sampleid", by.y = "cellLine", all.x = TRUE)
joined <- joined[, c("ctrp.sampleid", "accession", "AC")]

# compare each row of joined$accession and joined$AC
# return a data.table of the rows that are not equal
joined[accession != AC, .(ctrp.sampleid, accession, AC)]




