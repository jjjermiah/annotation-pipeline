library(readxl)
library(data.table)
library(Cellosaurus)


# get snakemake input and output files
inputfile <- snakemake@input[['dataset_annotation_file']]
outputfile <- snakemake@output[['sample_Cellosaurus_file']]

# get the wildcard value 
dataset_name <- snakemake@wildcards[['dataset']]

# GDSC_samples_excel <- "/home/bioinf/bhklab/jermiah/projects/annotationScripts/rawdata/GDSC/Cell_Lines_Details.xlsx"
excelsheets <- lapply(c(1, 2, 3), function(x){
    as.data.table(read_xlsx(path=inputfile, sheet=x))
})
names(excelsheets) <- c("cell_line_details", "COSMIC_tissue_classification", "Decode")


gdsc_sample_data <- excelsheets$cell_line_details 
gdsc_samples <- gdsc_sample_data$b
# remove the row where V1 == "TOTAL:"
gdsc_samples <- gdsc_samples[-which(gdsc_samples == "TOTAL:")]

# Use Cellosaurus functions from new package:  https://r.acidgenomics.com/packages/cellosaurus/articles/Cellosaurus.html
cellosaurus_object <- snakemake@input[['cellosaurus_object']]
object <- readRDS(cellosaurus_object)
# object <- Cellosaurus::Cellosaurus()


# get all cellosaurusdata into data.table
cellosaurusData <- as.data.table(object)

gdsc_mapped <- as.data.table(Cellosaurus::mapCells(object, cells = gdsc_samples),keep.rownames = T)
names(gdsc_mapped) <- c(paste0(dataset_name,".sampleid"), "cellosaurus.accession.ID") # set the names of the columns

any(is.na(gdsc_mapped$V2))
# Looks like they all map

# write ctrp_mapped to outputfile
fwrite(gdsc_mapped, outputfile, sep = "\t", quote = F, row.names = F, col.names = T)


