

# manually create the dictionary for each project's annotation file
# TODO:: moved into config.yaml file later
dataset_annotation_dict = {
    'ccle': {
        'sample': 'rawdata/ccle/Cell_lines_annotations_20181226.txt',
        'treatment': 'rawdata/ccle/CCLE_NP24.2009_Drug_data_2015.02.24.csv'
    },
    'ctrp': {
        'sample': 'rawdata/ctrp/v20.meta.per_cell_line.txt',
        'treatment': 'rawdata/ctrp/v20.meta.per_compound.txt'
    },
    'gdsc': {
        'sample': 'rawdata/gdsc/Cell_Lines_Details.xlsx',
        'treatment': 'rawdata/gdsc/screened_compounds_rel_8.4.csv'
    }
}
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


pubchemAnnotationQueries = [
    'ATC Code', 'NSC Number', 
    'ChEMBL ID', 'Drug Induced Liver Injury', 
    'CAS', 'FDA Approved Drugs', 
    'Names and Synonyms', 'Synonyms and Identifiers']

rule all:
    input:
        # sampleInfo_file = expand("results/{dataset}/sample_metadata.RDS", dataset = dataset_annotation_dict.keys()),
        pubchem_annotation_dbs = expand("metadata/pubchem_annotations_{query}.RDS", query = pubchemAnnotationQueries[2]),
        # treatmentInfo_file = expand("results/{dataset}/treatment_metadata.RDS", dataset = 'ctrp'),

rule getCellosaurus_annotation:
    input:
        dataset_annotation_file = lambda wildcards: dataset_annotation_dict[wildcards.dataset]['sample'],
        cellosaurus_object = "metadata/cellosaurus.RDS", 
    output:
        sample_Cellosaurus_file = "procdata/{dataset}/sample_Cellosaurus.csv",
    conda:
        "envs/sample_annotation.yaml",
    script:
        "scripts/getCellosaurus_{wildcards.dataset}.R"

rule getSampleMetadata:
    input:
        sample_Cellosaurus_file = "procdata/{dataset}/sample_Cellosaurus.csv",
        cellosaurus_object = "metadata/cellosaurus.RDS", 
        # sample_Cellosaurus_file = expand("procdata/{dataset}/sample_Cellosaurus.csv", dataset = dataset_annotation_dict.keys())
    output:
        sampleInfo_file = "results/{dataset}/sample_metadata.RDS",
    conda:
        "envs/sample_annotation.yaml",
    script:
        "scripts/CellosaurusAccessionID_to_metadata.R"

rule getCellosaurusObject:
    output:
        cellosaurus_object = "metadata/cellosaurus.RDS",
    conda:
        "envs/sample_annotation.yaml",
    script:
        "scripts/getCellosaurusObject.R"

rule getPubchemAnnotationQueryDatabase:
    output:
        pubchem_annotation_db = "metadata/pubchem_annotations_{query}.RDS",
    conda:
        "envs/sample_annotation.yaml",
    threads:
        8
    params:
        raw = "FALSE",
        rawAnnotationDT = "TRUE",
    script:
        "scripts/getPubchemAnnotations.R"

rule getPubchem_Annotation:
    input:
        dataset_annotation_file = lambda wildcards: dataset_annotation_dict[wildcards.dataset]['treatment']
    output:
        treatment_Pubchem_data = "results/{dataset}/pubchem_preproccessed.RDS",
    conda:
        "envs/sample_annotation.yaml",
    script:
        "scripts/getPubchem_{wildcards.dataset}.R"


# rule getTreatmentMetadata:
#     input:
#         treatment_Pubchem_data = "results/{dataset}/pubchem_preproccessed.RDS",
#     output:
#         treatmentInfo_file = "results/{dataset}/treatment_metadata.RDS",
#     conda:
#         "envs/sample_annotation.yaml",
#     script:
#         "scripts/getTreatmentMetadata_{wildcards.dataset}.R"
