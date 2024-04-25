## code to prepare `DATASET` dataset goes here

#wget https://www.science.org/doi/suppl/10.1126/science.adg7492/suppl_file/science.adg7492_data_s1_to_s9.zip
#unzip science.adg7492_data_s1_to_s9.zip
#cd science.adg7492_data_s1_to_s9.zip
#cp science.adg7492_data_s5.csv ClinVar_AM_supplement.csv

# Load in original CSV
clinvar <- read.csv("ClinVar_AM_supplement.csv")

# Separate UniProt ID and Protein Variant
clinvar <- tidyr::separate(clinvar, protein_variant, into = c("accession", "protein_variant"), sep = ":")

clinvar_data <- clinvar |>
    select(-('AlphaMissense')) |>
    rename(cv_variant_id = variant_id,
           cv_class = label,
           uniprot_id = accession) |>
    relocate('transcript_id', .after = 'uniprot_id')

usethis::use_data(clinvar_data, overwrite = TRUE)




