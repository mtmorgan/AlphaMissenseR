# Create GenomicRanges Objects from AlphaMissense Annotations

`to_GPos()` coerces a tibble derived from `am_data("hg19")` or
`am_data("hg38")` resources to GenomicRanges `GPos` objects.

## Usage

``` r
to_GPos(tbl)
```

## Arguments

- tbl:

  a tibble derived from `am_data("hg19")` or `am_data("hg38")`. The
  tibble must have columns `CHROM`, `POS`, and `genome`.

## Value

`to_GPos()` returns a `GPos` object, which can be used in the same was a
`GRanges` object for range-based filtering and annotation

## Examples

``` r
am_data("hg38") |>
    filter(CHROM == "chr2", POS < 10000000, REF == "G") |>
    select(-REF) |>
    to_GPos()
#> UnstitchedGPos object with 26751 positions and 6 metadata columns:
#>           seqnames       pos strand |         ALT  uniprot_id
#>              <Rle> <integer>  <Rle> | <character> <character>
#>       [1]     chr2     41613      * |           C      Q1W6H9
#>       [2]     chr2     41615      * |           C      Q1W6H9
#>       [3]     chr2     41615      * |           A      Q1W6H9
#>       [4]     chr2     41615      * |           T      Q1W6H9
#>       [5]     chr2     41618      * |           C      Q1W6H9
#>       ...      ...       ...    ... .         ...         ...
#>   [26747]     chr2   9999028      * |           A      Q9NZI5
#>   [26748]     chr2   9999028      * |           T      Q9NZI5
#>   [26749]     chr2   9999029      * |           C      Q9NZI5
#>   [26750]     chr2   9999029      * |           A      Q9NZI5
#>   [26751]     chr2   9999029      * |           T      Q9NZI5
#>                transcript_id protein_variant am_pathogenicity          am_class
#>                  <character>     <character>        <numeric>       <character>
#>       [1]  ENST00000327669.5           R321G           0.0848     likely_benign
#>       [2]  ENST00000327669.5           S320C           0.0946     likely_benign
#>       [3]  ENST00000327669.5           S320F           0.1593     likely_benign
#>       [4]  ENST00000327669.5           S320Y           0.1429     likely_benign
#>       [5]  ENST00000327669.5           P319R           0.0765     likely_benign
#>       ...                ...             ...              ...               ...
#>   [26747] ENST00000324907.14           G581R           0.9997 likely_pathogenic
#>   [26748] ENST00000324907.14           G581W           0.9996 likely_pathogenic
#>   [26749] ENST00000324907.14           G581A           0.9970 likely_pathogenic
#>   [26750] ENST00000324907.14           G581E           0.9995 likely_pathogenic
#>   [26751] ENST00000324907.14           G581V           0.9994 likely_pathogenic
#>   -------
#>   seqinfo: 25 sequences (1 circular) from hg38 genome
```
