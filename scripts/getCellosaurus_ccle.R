library(data.table)
library(Cellosaurus)

path <- "/home/bioinf/bhklab/jermiah/projects/annotationScripts/rawdata/ccle/Cell_lines_annotations_20181226.txt"
dataset_name <- "ccle"
# get snakemake input and output files
inputfile <- snakemake@input[['dataset_annotation_file']]
outputfile <- snakemake@output[['sample_Cellosaurus_file']]

# get the wildcard value 
dataset_name <- snakemake@wildcards[['dataset']]

ccle <- fread(path, header = T, sep = "\t")

# Use Cellosaurus functions from new package:  https://r.acidgenomics.com/packages/cellosaurus/articles/Cellosaurus.html
cellosaurus_object <- snakemake@input[['cellosaurus_object']]
object <- readRDS(cellosaurus_object)
# object <- Cellosaurus::Cellosaurus()


# get all cellosaurusdata into data.table
cellosaurusData <- as.data.table(object)

# First use depmapID to map to cellosaurus.accession.ID
mapped_depmap <- as.data.table(mapCells(object, cells = ccle[!is.na(ccle$depMapID),depMapID]), keep.rownames=T)
# merge mapped_depmap with ccle to get the CCLE_ID 
mapped_depmap <- merge(mapped_depmap, ccle, by.x = "V1", by.y = "depMapID", all.y=T)[,c("CCLE_ID", "V2")]

# get CCLE_IDs that didn't map
missing.CCLE_ID <- mapped_depmap[is.na(V2), c("CCLE_ID")]
missing.CCLE_ID[,CCLE_sample_name := tstrsplit(x=CCLE_ID, split="_")[1][1],]
mapped_missing <- as.data.table(mapCells(object, cells = missing.CCLE_ID$CCLE_sample_name), keep.rownames=T)
mapped_missing <- merge(missing.CCLE_ID, mapped_missing, by.x = "CCLE_sample_name", by.y = "V1", all.y=T)[,c("CCLE_ID", "V2")]

ccle_mapped <- rbindlist(list(mapped_depmap[!is.na(V2)], mapped_missing[!is.na(V2)]))

# # testing using Name
# mapped_name <- as.data.table(mapCells(object, cells = ccle[!is.na(Name), Name]), keep.rownames=T)

names(ccle_mapped) <- c(paste0(dataset_name,".sampleid"), "cellosaurus.accession.ID") # set the names of the columns

fwrite(ccle_mapped, outputfile, sep = "\t", quote = F, row.names = F, col.names = T)