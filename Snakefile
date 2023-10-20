
from os.path import join
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

pubchemAnnotationQueries = ['ATC Code', 'NSC Number', 'ChEMBL ID', 'Drug Induced Liver Injury', 'CAS']
DATASETS = "ctrp" #["ccle", "ctrp", "gdsc"]
rule all:
    input:
        treatmentMetadata = expand("procdata/{dataset}/pubchem_annotations_all.csv", dataset = DATASETS), #"results/ctrp/treatment_metadata.RDS",
        sampleMetadata = expand("results/ctrp/sample_metadata.RDS", dataset = "ctrp"),

rule combineSampleCellosaurusMetadata:
    input:
        sample_Cellosaurus_file = "procdata/{dataset}/sample_Cellosaurus.csv",
        cellosaurus_object = "metadata/cellosaurus.RDS", 
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
        "scripts/getCellosaurus/getCellosaurusObject.R"

rule mapSampleNamesToCellosaurusAccessionID:
    input:
        dataset_annotation_file = lambda wildcards: dataset_annotation_dict[wildcards.dataset]['sample'],
        cellosaurus_object = "metadata/cellosaurus.RDS", 
    output:
        sample_Cellosaurus_file = "procdata/{dataset}/sample_Cellosaurus.csv",
    conda:
        "envs/sample_annotation.yaml",
    script:
        "scripts/getCellosaurus/getCellosaurus_{wildcards.dataset}.R"

rule combinePubchemAnnotations:
    input:
        pubchem_annotation_db = expand("procdata/{dataset}/pubchem_annotations_{annotationType}.RDS", dataset = DATASETS,  annotationType = pubchemAnnotationQueries),
        pubchem_properties = "procdata/{dataset}/pubchem_properties.RDS"
    output:
        pubchemTreatmentAll = "procdata/{dataset}/pubchem_annotations_all.csv"
    conda:
        "envs/sample_annotation.yaml",
    script:
        "scripts/getPubChem/combinePubchemAnnotations.R"

rule getExternalAnnotationFromCID:
    input:
        treatment_Pubchem_data = "procdata/{dataset}/pubchem_preproccessed.RDS",
    output:
        pubchem_annotation_dbs = "procdata/{dataset}/pubchem_annotations_{annotationType}.RDS"
    threads: 10
    conda:
        "envs/sample_annotation.yaml",
    script:
        "scripts/getPubChem/getPubchemAnnotation.R"

rule getPropertiesFromCID:
    input:
        treatment_Pubchem_data = "procdata/{dataset}/pubchem_preproccessed.RDS",
    output:
        pubchem_properties = "procdata/{dataset}/pubchem_properties.RDS"
    threads: 10
    conda:
        "envs/sample_annotation.yaml",
    script:
        "scripts/getPubChem/getPubchemProperties.R"

rule mapTreatmentNamesToPubchemCID:
    input:
        dataset_annotation_file = lambda wildcards: dataset_annotation_dict[wildcards.dataset]['treatment']
    output:
        treatment_Pubchem_data = "procdata/{dataset}/pubchem_preproccessed.RDS",
    conda:
        "envs/sample_annotation.yaml",
    script:
        "scripts/getPubChem/getPubChem_{wildcards.dataset}.R"

rule getChEMBLAnnotations:
    input:
        CIDtoCHEMBL = expand("procdata/{dataset}/pubchem_annotations_{annotationType}.RDS", dataset = DATASETS, annotationType = ['ChEMBL ID']),
    output:
        allChEMBLMechanisms = "procdata/{dataset}/ChEMBL_mechanisms.csv"
    conda:
        "envs/sample_annotation.yaml",
    script:
        "scripts/getChEMBL/getChEMBLAnnotations.R"

# rule getPubchemAnnotationQueryDatabase:
#     output:
#         pubchem_annotation_db = expand("procdata/pubchem_annotations_{query}.RDS", query = pubchemAnnotationQueries)
#     conda:
#         "envs/sample_annotation.yaml",
#     threads:
#         10
#     params:
#         raw = "FALSE",
#         rawAnnotationDT = "TRUE",
#         maxPages = "5",
#         verbose = "TRUE"
#     script:
#         "scripts/getPubchemAnnotations.R"