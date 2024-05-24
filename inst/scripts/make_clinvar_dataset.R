#wget https://www.science.org/doi/suppl/10.1126/science.adg7492/suppl_file/science.adg7492_data_s1_to_s9.zip
#unzip science.adg7492_data_s1_to_s9.zip
#cd science.adg7492_data_s1_to_s9
#cp science.adg7492_data_s5.csv .
#gzip science.adg7492_data_s5.csv
library(tidyr)
library(dplyr)

# Load in original CSV
fpath <- system.file(
    "extdata",
    "science.adg7492_data_s5.csv.gz",
    package = "AlphaMissenseR"
)

clinvar <- read.csv(fpath)

# Separate UniProt ID and Protein Variant
clinvar <-
    clinvar |>
    tidyr::separate(
        .data$protein_variant,
        into = c("uniprot_id", "protein_variant"),
        sep = ":"
    )

clinvar_data <-
    clinvar |>
    select(-('AlphaMissense')) |>
    rename(cv_variant_id = variant_id,
           cv_class = label) |>
    relocate('transcript_id', .after = 'uniprot_id')

usethis::use_data(clinvar_data, overwrite = TRUE)
