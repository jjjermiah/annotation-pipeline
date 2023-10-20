## ------------------- Parse Snakemake Object ------------------- ##
if(exists("snakemake")){
    inputFile <- snakemake@input[['CIDtoCHEMBL']]
    inputPubchemProperties <- snakemake@input[['pubchem_properties']]
    annotationTypes <- sapply(inputFiles, function(x){
        file <- paste((basename(x)), collapse=" ")      # remove the directories from path
        file <- strsplit(file, "_")[[1]]                # split on underscore 
        strsplit(file[[length(file)]], "\\.")[[1]][[1]] # get last element and remove the extension
    })
    outputFile <- snakemake@output[['pubchemTreatmentAll']]
    save.image()
}

