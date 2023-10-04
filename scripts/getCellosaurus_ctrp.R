#  GetCellosaurus_ctrp.R 

library(data.table)
library(Cellosaurus)

# inputfile <- "/home/bioinf/bhklab/jermiah/projects/annotationScripts/rawdata/CTRPv2/v20.meta.per_cell_line.txt"

# get snakemake input and output files
inputfile <- snakemake@input[['dataset_annotation_file']]
outputfile <- snakemake@output[['sample_Cellosaurus_file']]

# get the wildcard value 
dataset_name <- snakemake@wildcards[['dataset']]

########
# LOAD FILES
ctrpcells <- fread(inputfile, header = T, sep = "\t")

# Use Cellosaurus functions from new package:  https://r.acidgenomics.com/packages/cellosaurus/articles/Cellosaurus.html
cellosaurus_object <- snakemake@input[['cellosaurus_object']]
object <- readRDS(cellosaurus_object)
# object <- Cellosaurus::Cellosaurus()


# get all cellosaurusdata into data.table
cellosaurusData <- as.data.table(object)

# use mapping function, resulting ctrp_mapped is a named 'character' converted to a data.table

ctrp_mapped <- as.data.table(Cellosaurus::mapCells(object, cells = ctrpcells$ccl_name), keep.rownames = T)
names(ctrp_mapped) <- c(paste0(dataset_name,".sampleid"), "cellosaurus.accession.ID") # set the names of the columns

# write ctrp_mapped to outputfile
fwrite(ctrp_mapped, outputfile, sep = "\t", quote = F, row.names = F, col.names = T)