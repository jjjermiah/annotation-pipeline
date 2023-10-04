library(data.table)
library(Cellosaurus)
# library(tidyverse)


# install R-package from github
# repo: bhklab/AnnotationGx
# branch: "jermiah"
# commitid: 634db68
#
# install.packages("devtools")
# devtools::install_github("bhklab/AnnotationGx@634db68", ref = "jermiah")
devtools::install_github("bhklab/AnnotationGx@634db68", ref = "jermiah", dependencies = TRUE, upgrade="never")
library(AnnotationGx)

AnnotationGx::getCellosaurusAPI
# get snakemake input and output files
inputfiles <- snakemake@input[['sample_Cellosaurus_file']]
cellosaurus_object <- snakemake@input[['cellosaurus_object']]
outputfile <- snakemake@output[['sampleInfo_file']]

# get the wildcard value 
dataset_name <- snakemake@wildcards[['dataset']]


# read in the sample cellosaurus file
sample_cellosaurus <- fread(inputfiles[1], sep = "\t", header = T, stringsAsFactors = F)
# print(sample_cellosaurus)


# sample_cellosaurus <-fread("/home/bioinf/bhklab/jermiah/projects/annotationScripts/procdata/ctrp/sample_Cellosaurus.csv", sep = "\t", header = T, stringsAsFactors = F)
object <- readRDS(cellosaurus_object)
# get all cellosaurusdata into data.table
cellosaurusData <- as.data.table(object)

columns <- c(
    "cellLineName", "synonyms", "accession", "misspellings",
    # "diseases", 
    "category", "samplingSite", "isCancer", 
    "oncotreeName", "oncotreeTissue", "oncotreeLevel",
    "depmapId", "atccId", "ncbiTaxonomyId", "sangerModelId", "secondaryAccession",
    "sexOfCell", "ageAtSampling"
)

# subset cellosaurusData to only the columns we want
cellosaurusData <- cellosaurusData[, columns, with = F]

merged_dt <- merge(sample_cellosaurus, cellosaurusData, by.x = "cellosaurus.accession.ID", by.y = "accession", all.x = TRUE)

agx_results <- 
    AnnotationGx::getCellosaurusAPI(
        cl_names = unique(merged_dt[cellosaurus.accession.ID!="", cellosaurus.accession.ID]), 
        fields = c("din"),
        GETfxn = "search/cell-line?",
        querydomain = "ac:"
    )
cleaned <- AnnotationGx::cleanCellosaurusResponse(agx_results,GETfxn = "search/cell-line?")
cleaned[, c("NCIt.ID", "NCIt.Name"):= tstrsplit(DI, split = ";")[2:3]][, DI:=NULL]

merged_dt <-  merge(merged_dt, cleaned,  by.x="cellosaurus.accession.ID", by.y="cellLine", all.x=T, sort=T)

# write merged_dt  to outputfile
# convert any empty or NA cells to a character string of "NA"
merged_dt[is.na(merged_dt)] <- "NA"
saveRDS(merged_dt, outputfile)

# results <- as.data.table(
#     tidyr::unnest_wider(
#         tidyr::unnest_wider(
#             merged_dt, 
#             diseases
#         ),
#         col = "NCIt",
#         names_sep = "_"
#     )
# )


##############################################################################################
# Get all cross references

# & get all the Sources from the resulting column names excluding "accession"
# crossReference_dt <- as.data.table(
#     tidyr::unnest_wider(
#         cellosaurusData[,c("accession", "crossReferences")],
#         crossReferences
#     )
# )
# crossReferenceSources <- colnames(crossReference_dt)[-1]


##############################################################################################



