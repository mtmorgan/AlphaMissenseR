# AlphaMissense Data for *R* / *Bioconductor*

The AlphaMissense
[publication](https://www.science.org/doi/epdf/10.1126/science.adg7492)
outlines how a variant of AlphaFold / DeepMind was used to predict
missense variant pathogenicity. Supporting data on
[Zenodo](https://zenodo.org//record/8360242) include, for instance 70+M
variants across hg19 and hg38 genome builds. The AlphaMissense package
allows ready access to the data, downloading individual files to DuckDB
databases for ready exploration and integration into *R* and
*Bioconductor* worksflows.

## Installation

Install the package from Bioconductor or GitHub, ensuring correct
*Bioconductor* dependencies.

When the package is available on *Bioconductor*, use

``` r
if (!"BiocManager" %in% rownames(installed.packages()))
    install.packages("BiocManager", repos = "https://cloud.R-project.org")
if (BiocManager::version() >= "3.19") {
    BiocManager::install("AlphaMissenseR")
} else {
    stop(
        "'AlphaMissenseR' requires Bioconductor version 3.19 or later, ",
        "install from GitHub?"
    )
}
```

Use the pre-release or devel version with

``` r
if (!"remotes" %in% rownames(installed.packages()))
    install.packages("remotes", repos = "https://cloud.R-project.org")
remotes::install_github(
    "mtmorgan/AlphaMissenseR",
    repos = BiocManager::repositories()
)
```

Load the library.

``` r
library(AlphaMissenseR)
```

## Next steps

- Visit the
  [Introduction](https://mtmorgan.github.io/AlphaMissenseR/articles/introduction.md)
  to learn more about accessing AlphaMissense data in *R*.
- The [AlphaFold
  Integration](https://mtmorgan.github.io/AlphaMissenseR/articles/alphafold.md)
  article shows how missense effects can be plotted on AlphaFold (or
  other) protein structures.
- Use [ClinVar
  Integration](https://mtmorgan.github.io/AlphaMissenseR/articles/clinvar.md)
  to compare AlphaMissense and ClinVar classifications.
- See [Benchmarking with
  ProteinGym](https://mtmorgan.github.io/AlphaMissenseR/articles/benchmarking.md)
  for assessing AlphaMissense predictions.
- Troubleshoot common problems with [Issues &
  Solutions](https://mtmorgan.github.io/AlphaMissenseR/articles/issues.md).
