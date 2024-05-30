#lynx -accept_all_cookies=TRUE https://www.science.org/doi/suppl/10.1126/science.adg7492/suppl_file/science.adg7492_data_s1_to_s9.zip

# Unzip directory of supplmental files to access separate datasets
unzip('inst/extdata/science.adg7492_data_s1_to_s9.zip', exdir = 'inst/extdata/')

# Compress the ClinVar dataset size in order to be included in R package
R.utils::gzip('inst/extdata/science.adg7492_data_s5.csv')

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
    mutate(cv_class = as.factor(cv_class)) |>
    relocate('transcript_id', .after = 'uniprot_id')
