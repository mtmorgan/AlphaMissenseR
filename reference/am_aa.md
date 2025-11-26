# Common Amino Acid-level Transformations

`am_aa_pos()` separates `protein_variant` columns into amino acide
'pos', 'ref', and 'alt' columns.

`am_aa_pathogenicity()` summarizes pathogenicity scores at each protein
amino acid position.

## Usage

``` r
am_aa_pos(tbl)

am_aa_pathogenicity(tbl)
```

## Arguments

- tbl:

  a tibble, usually derived from `am_data("aa_substitutions")`,
  \`am_data("hg38"), etc. See details.

## Value

`am_aa_pos()` returns the original table with additional columns

- `aa_pos`: the position of the protein variant, as an
  [`integer()`](https://rdrr.io/r/base/integer.html).

- `aa_ref`: the single-character reference amino acid in the protein
  variant.

- `aa_alt`: the single-character alternate amino acid in the protein
  variant.

`am_aa_pathogenicity()` returns a tibble with columns

- `uniprot_id`, `aa_pos`, `aa_ref`: the UniProt id, and the position and
  reference amino acid being summarized

- `aa_pathogenicity_n`, `aa_pathogenicity_mean`,
  `aa_pathogenicity_median`, `aa_pathogenicity_min`,
  `aa_pathogenicity_max`: the number, average, median, minimum, and
  maximum of the pathogenicity scores at each amino acid position.

- `aa_pathogenicity_mode`: the modal `am_class` at the amino acid
  position, as a factor. Tied mode is assigned to lower pathogenicity.

## Details

`tbl` is
[`collect()`](https://dplyr.tidyverse.org/reference/compute.html)ed
before computation, so all rows must fit into memory.

For `am_aa_pos()`, `tbl` must contain a column `protein_variant` with
entries in the form `"Q465H"`, as in the AlphaMissense data.

For `am_aa_pathogenicity()`, `tbl` must contain columns `uniprot_id`,
`protein_variant`, `am_pathogenicity` and `am_class`. If `am_pos` and
friends are not already calculated, then `am_aa_pos()` is called.

## Examples

``` r
P35557 <-
    am_data("hg38") |>
    filter(uniprot_id %in% "P35557")

am_aa_pos(P35557)
#> # A tibble: 3,085 × 13
#>    CHROM      POS REF   ALT   genome uniprot_id transcript_id    protein_variant
#>    <chr>    <dbl> <chr> <chr> <chr>  <chr>      <chr>            <chr>          
#>  1 chr7  44145139 C     A     hg38   P35557     ENST00000403799… Q465H          
#>  2 chr7  44145139 C     G     hg38   P35557     ENST00000403799… Q465H          
#>  3 chr7  44145140 T     C     hg38   P35557     ENST00000403799… Q465R          
#>  4 chr7  44145140 T     A     hg38   P35557     ENST00000403799… Q465L          
#>  5 chr7  44145140 T     G     hg38   P35557     ENST00000403799… Q465P          
#>  6 chr7  44145141 G     C     hg38   P35557     ENST00000403799… Q465E          
#>  7 chr7  44145141 G     T     hg38   P35557     ENST00000403799… Q465K          
#>  8 chr7  44145143 C     G     hg38   P35557     ENST00000403799… G464A          
#>  9 chr7  44145143 C     T     hg38   P35557     ENST00000403799… G464D          
#> 10 chr7  44145143 C     A     hg38   P35557     ENST00000403799… G464V          
#> # ℹ 3,075 more rows
#> # ℹ 5 more variables: am_pathogenicity <dbl>, am_class <chr>, aa_pos <int>,
#> #   aa_ref <chr>, aa_alt <chr>

am_aa_pos(P35557) |>
    select(
        uniprot_id, POS, REF, ALT, protein_variant,
        starts_with("aa_"), am_pathogenicity, am_class
    ) |>
    arrange(aa_pos)
#> # A tibble: 3,085 × 10
#>    uniprot_id      POS REF   ALT   protein_variant aa_pos aa_ref aa_alt
#>    <chr>         <dbl> <chr> <chr> <chr>            <int> <chr>  <chr> 
#>  1 P35557     44188949 A     C     L2R                  2 L      R     
#>  2 P35557     44188949 A     T     L2Q                  2 L      Q     
#>  3 P35557     44188949 A     G     L2P                  2 L      P     
#>  4 P35557     44188950 G     T     L2M                  2 L      M     
#>  5 P35557     44188950 G     C     L2V                  2 L      V     
#>  6 P35557     44188945 G     C     D3E                  3 D      E     
#>  7 P35557     44188945 G     T     D3E                  3 D      E     
#>  8 P35557     44188946 T     G     D3A                  3 D      A     
#>  9 P35557     44188946 T     C     D3G                  3 D      G     
#> 10 P35557     44188946 T     A     D3V                  3 D      V     
#> # ℹ 3,075 more rows
#> # ℹ 2 more variables: am_pathogenicity <dbl>, am_class <chr>


am_aa_pathogenicity(P35557)
#> # A tibble: 464 × 9
#>    uniprot_id aa_pos aa_ref aa_pathogenicity_n aa_pathogenicity_mean
#>    <chr>       <int> <chr>               <int>                 <dbl>
#>  1 P35557          2 L                       5                0.0818
#>  2 P35557          3 D                       8                0.184 
#>  3 P35557          4 D                       8                0.147 
#>  4 P35557          5 R                       6                0.250 
#>  5 P35557          6 A                       6                0.138 
#>  6 P35557          7 R                       7                0.257 
#>  7 P35557          8 M                       9                0.142 
#>  8 P35557          9 E                       7                0.212 
#>  9 P35557         10 A                       6                0.142 
#> 10 P35557         11 A                       6                0.142 
#> # ℹ 454 more rows
#> # ℹ 4 more variables: aa_pathogenicity_median <dbl>,
#> #   aa_pathogenicity_min <dbl>, aa_pathogenicity_max <dbl>,
#> #   aa_pathogenicity_mode <fct>
```
