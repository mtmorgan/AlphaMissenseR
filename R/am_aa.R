#' @rdname am_aa
#'
#' @title Common Amino Acid-level Transformations
#'
#' @description `am_aa_pos()` separates `protein_variant` columns into
#'     amino acide 'pos', 'ref', and 'alt' columns.
#'
#' @param tbl a tibble, usually derived from
#'     `am_data("aa_substitutions")`, `am_data("hg38"), etc. See
#'     details.
#'
#' @details
#'
#' `tbl` is `collect()`ed before computation, so all rows must fit
#' into memory.
#'
#' For `am_aa_pos()`, `tbl` must contain a column `protein_variant`
#' with entries in the form `"Q465H"`, as in the AlphaMissense data.
#'
#' @return
#'
#' `am_aa_pos()` returns the original table with additional columns
#'
#' - `aa_pos`: the position of the protein variant,
#'   as an `integer()`.
#' - `aa_ref`: the single-character reference amino acid in the
#'   protein variant.
#' - `aa_alt`: the single-character alternate amino acid in the
#'   protein variant.
#'
#' @examples
#'
#' P35557 <-
#'     am_data("hg38") |>
#'     filter(uniprot_id %in% "P35557")
#'
#' am_aa_pos(P35557)
#'
#' am_aa_pos(P35557) |>
#'     select(
#'         uniprot_id, POS, REF, ALT, protein_variant,
#'         starts_with("aa_"), am_pathogenicity, am_class
#'     ) |>
#'     arrange(aa_pos)
#'
#' @export
am_aa_pos <-
    function(tbl)
{
    required_columns <- "protein_variant"
    stopifnot(
        inherits(tbl, "tbl"),
        all(required_columns %in% colnames(tbl))
    )

    pattern <- "([A-Z])([[:digit:]]+)([A-Z])"
    result <-
        tbl |>
        collect() |>
        mutate(
            aa_pos = as.integer(sub(pattern, "\\2", .data$protein_variant)),
            aa_ref = sub(pattern, "\\1", .data$protein_variant),
            aa_alt = sub(pattern, "\\3", .data$protein_variant)
        )

    result
}

am_aa_class_mode <-
    function(x)
{
    am_class_levels <- c("likely_benign", "ambiguous", "likely_pathogenic")
    if (length(x)) {
        n <- tabulate(match(x, am_class_levels))
        ## ties go to lower-pathogenicity level
        x <- am_class_levels[which.max(n)]
    }
    ## return as factor
    factor(x, am_class_levels)
}

#' @rdname am_aa
#'
#' @description `am_aa_pathogenicity()` summarizes pathogenicity
#'     scores at each protein amino acid position.
#'
#' @details
#'
#' For `am_aa_pathogenicity()`, `tbl` must contain columns
#' `uniprot_id`, `protein_variant`, `am_pathogenicity` and
#' `am_class`. If `am_pos` and friends are not already calculated,
#' then `am_aa_pos()` is called.
#'
#' @return
#'
#' `am_aa_pathogenicity()` returns a tibble with columns
#'
#' - `uniprot_id`, `aa_pos`, `aa_ref`: the UniProt id, and the
#'   position and reference amino acid being summarized
#' - `aa_pathogenicity_n`, `aa_pathogenicity_mean`,
#'   `aa_pathogenicity_median`, `aa_pathogenicity_min`,
#'   `aa_pathogenicity_max`: the number, average, median, minimum, and
#'   maximum of the pathogenicity scores at each amino acid position.
#' - `aa_pathogenicity_mode`: the modal `am_class` at the amino acid
#'   position, as a factor. Tied mode is assigned to lower
#'   pathogenicity.
#'
#' @importFrom dplyr group_by summarize ungroup n
#'
#' @importFrom stats median
#'
#' @examples
#'
#' am_aa_pathogenicity(P35557)
#'
#' @export
am_aa_pathogenicity <-
    function(tbl)
{
    required_columns <- c(
        "uniprot_id", "protein_variant", "am_pathogenicity", "am_class"
    )
    stopifnot(
        inherits(tbl, "tbl"),
        all(required_columns %in% colnames(tbl))
    )

    if (!"aa_pos" %in% colnames(tbl))
        tbl <- am_aa_pos(tbl)

    aa_summary <-
        tbl |>
        group_by(.data$uniprot_id, .data$aa_pos, .data$aa_ref) |>
        summarize(
            aa_pathogenicity_n = n(),
            aa_pathogenicity_mean = mean(.data$am_pathogenicity),
            aa_pathogenicity_median = median(.data$am_pathogenicity),
            aa_pathogenicity_min = min(.data$am_pathogenicity),
            aa_pathogenicity_max = max(.data$am_pathogenicity),
            aa_pathogenicity_mode = am_aa_class_mode(.data$am_class)
        ) |>
        ungroup()

    aa_summary
}
