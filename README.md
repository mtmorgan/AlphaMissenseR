
# AlphaMissense for *R* / *Bioconductor*

<!-- badges: start -->
<!-- badges: end -->

The AlphaMissense [publication][Science] outlines how a variant of
AlphaFold / DeepMind was used to predict missense variant
pathogenicity. Supporting data on [Zenodo][] include, for instance
70+M variants across hg19 and hg38 genome builds. The AlphaMissense
package allows ready access to the data, downloading individual files
to DuckDB databases for ready exploration and integration into *R* and
*Bioconductor* worksflows.

[Science]: https://www.science.org/doi/epdf/10.1126/science.adg7492
[Zenodo]: https://zenodo.org//record/8360242

## Installation

Install the package from GitHub, ensuring correct *Bioconductor*
dependencies.

``` r
if (!"BiocManager" %in% rownames(installed.packages()))
    install.packages("BiocManager", repos = "https://cran.r-project.org")

remotes::install_github(
    "mtmorgan/AlphaMissense",
    repos = BiocManager::repositories()
)
```

Load the library.

``` r
library(AlphaMissense)
```

## Next steps

Visit the [introductory article][intro] to learn more about using this
package.

[intro]: https://mtmorgan.github.io/AlphaMissense/articles/introduction.html
