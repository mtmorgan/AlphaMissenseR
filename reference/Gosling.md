# Plot a GPos or GRanges Object Using Shiny and Gosling

`plot_granges()` creates a Shiny app that displays a Gosling plot of a
given GRanges object. It visualizes genomic ranges with both rectangle
and point representations, and allows for customization of the plot
title and subtitle.

## Usage

``` r
gosling_plot(
  granges,
  title = "GRanges Plot",
  subtitle = "Stacked nucleotide example",
  plot_type = c("bar", "lollipop"),
  create_app = interactive()
)
```

## Arguments

- granges:

  A
  [`GenomicRanges::GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  (e.g.,
  [`GenomicRanges::GPos`](https://rdrr.io/pkg/GenomicRanges/man/GPos-class.html))
  object containing the genomic ranges to be plotted.

- title:

  character(1) Title of the plot. Default is "GRanges Plot".

- subtitle:

  character(1) The subtitle of the plot. Default is "Stacked nucleotide
  example".

- plot_type:

  character(1) Select the type of gosling plot. Default is "bar";
  available types are described in Details.

- create_app:

  logical(1) Produce a Shiny app for Gosling visualization. Default
  `TRUE` when used interactively.

## Value

`gosling_plot()` with `create_app = TRUE` returns a `shinyApp` object
that, when run, displays the Gosling plot. When `create_app = FALSE`,
the (invisible) return value represents the gosling object to be used in
Shiny apps.

## Details

The function supports 2 types of plots as selected through the
`plot_type` argument: (1) `"bar"`: a barplot-like view of the
pathogenicity score at each position, similar to a sequencing coverage
plot, and (2) `"lollipop"`: a lollipop plot focusing on the
pathogenicity classification (ambiguous, benign, pathogenic) at each
position. It requires that the
[`GenomicRanges::GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
object has the following metadata columns: 'am_class' (effect
classification),'am_pathogenicity' (pathogenicity score), 'ALT'
(alternative allele) and 'REF' (reference allele).

## Note

This function requires the `shiny`, `shiny.gosling`, and `GenomicRanges`
packages to be installed.

## Examples

``` r
## Create a sample GRanges object from AlphamissenseR
gpos <-
    am_data("hg38") |>
    filter(uniprot_id == "Q1W6H9") |>
    to_GPos()

## Plot `gpos` using shiny.gosling
if (requireNamespace("shiny.gosling", quietly = TRUE)) {
   gosling_plot(
       gpos, title = "Q1W6H_track", subtitle = "bar plot example",
       plot_type = "bar", create_app = FALSE
   )
}
```
