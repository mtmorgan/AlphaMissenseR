# Download Supplemental Materials from AlphaMissense publication
# browseURL('https://www.science.org/doi/suppl/10.1126/science.adg7492/suppl_file/science.adg7492_data_s1_to_s9.zip'

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
