
# AlphaMissense Data for *R* / *Bioconductor*

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

Install the package from Bioconductor or GitHub, ensuring correct
*Bioconductor* dependencies.

When the package is available on *Bioconductor*, use

```{r install-Bioconductor, eval = FALSE}
BiocManager::install("AlphaMissenseData")
```

Use the pre-release or devel version with

```{r install-devel, eval = FALSE}
remotes::install_github(
    "mtmorgan/AlphaMissenseData",
    repos = BiocManager::repositories()
)
```

Load the library.

``` r
library(AlphaMissenseData)
```

## Next steps

Visit the [introductory article][intro] to learn more about using this
package.

[intro]: https://mtmorgan.github.io/AlphaMissenseData/articles/introduction.html
