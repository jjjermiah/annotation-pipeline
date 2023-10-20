library(Cellosaurus)

outputfile <- snakemake@output[['cellosaurus_object']]
object <- Cellosaurus::Cellosaurus()

# write object to outputfile as .RDS file
saveRDS(object, outputfile)
