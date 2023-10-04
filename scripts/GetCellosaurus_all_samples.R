library(data.table)
library(Cellosaurus)


# get snakemake input and output files
inputfile <- snakemake@input[['dataset_annotation_file']]
outputfile <- snakemake@output[['sample_Cellosaurus_file']]

# sample_github_csv <- "/home/bioinf/bhklab/jermiah/projects/annotationScripts/rawdata/cell_annotation_all.csv"
sample_github_csv <- as.data.table(read.csv(sample_github_csv, header = T, sep = ","))

# Use Cellosaurus functions from new package:  https://r.acidgenomics.com/packages/cellosaurus/articles/Cellosaurus.html
object <- Cellosaurus::Cellosaurus()

# get all cellosaurusdata into data.table
cellosaurusData <- as.data.table(object)

all_mapped_cells <- as.data.table(Cellosaurus::mapCells(object, cells = sample_github_csv$unique.cellid),keep.rownames = T)
names(all_mapped_cells) <- c("unique.sampleid", "cellosaurus.accession.ID") # set the names of the columns

merged_dt <- merge(cellosaurusData, mapped_cells, by.x = "accession", by.y = "cellosaurus.accession.ID", all.y = TRUE)

