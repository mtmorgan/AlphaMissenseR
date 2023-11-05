ALPHAFOLD_PREDICTION <-
    "https://alphafold.ebi.ac.uk/api/prediction/{{qualifier}}"

af_prediction <-
    function(accession)
{
    nmspcs <-
        requireNamespace("httr", quietly = TRUE) &&
        requireNamespace("tidyr", quietly = TRUE)
    if (!nmspcs)
        stop("`af_preditions()` requires packages 'httr', 'tidyr'")
    url <- whisker.render(ALPHAFOLD_PREDICTION, c(qualifier = accession))
    response <- httr::GET(url)
    if (httr::status_code(response) < 400L) {
        record <- httr::content(response)
        tibble(record) |>
            tidyr::unnest_wider(.data$record)
    } else {
        spdl::debug(
            "af_prediction() accession '{}' status_code '{}'",
            accession, httr::status_code(response)
        )
        NA
    }
}

#' @rdname AlphaFold
#'
#' @title AlphaFold Protein Structure Retrieval and Display
#'
#' @description `af_predictions()` retrieves information about
#'     AlphaFold predictions associated with UniProt accession
#'     identifiers.
#'
#' @param uniprot_ids `character()` UniProt accession identifiers
#'     (`uniprot_id` in AlphaMissense tables).
#'
#' @details
#'
#' `af_predictions()` queries the `prediction` endpoint of the
#' AlphaFold API described at <https://alphafold.ebi.ac.uk/api-docs>.
#'
#' @return
#'
#' `af_predictions()` returns a tibble. Each row represents the
#' AlphaFold prediction associated with the corresponding uniprot
#' accession. Columns include:
#'
#' - `entryId`: AlphaFold identifier.
#' - `gene`: gene symbol corresponding to UniProt protein.
#' - `uniprotAccession`, `uniprotId`, `uniprotDescription`:
#'   UniProt characterization. AlphaMissense's `uniprot_id` is
#'   AlphaFold's `uniprotAccession`.
#' - `taxId`, `organismScientificName`: Organism information.
#' - `uniprotStart`, `uniprotEnd`, `uniprotSequence`:
#'   protein sequence information.
#' - `modelCreatedDate`, `latestVersion`, `allVersions`,
#'   `isReviewed`, `isReferenceProteome`: AlphaFold provenance
#'   information.
#' - `cifUrl`, `bcifUrl`, `pdbUrl`:
#'   URLs to AlphaFold 3-dimensional molecular representations.
#' - `paeImageUrl`, `paeDocUrl`: 'Predicted Aligned Error' heat map
#'   and underlying data. These can be used to assess the confidence
#'   in relative orientation of residues in different domains, as
#'   described in part in the AlphaFold FAQ
#'   <https://alphafold.ebi.ac.uk/faq>
#'
#' @examples
#'
#' ## af_predictions
#'
#' uniprot_ids <-
#'     am_data("aa_substitutions") |>
#'     dplyr::filter(uniprot_id %like% "P3555%") |>
#'     dplyr::distinct(uniprot_id) |>
#'     pull(uniprot_id)
#' af_predictions(uniprot_ids)
#'
#' @importFrom whisker whisker.render
#'
#' @importFrom dplyr bind_rows
#'
#' @importFrom utils head
#'
#' @export
af_predictions <-
    function(uniprot_ids)
{
    stopifnot(
        length(uniprot_ids) >= 1L,
        is_character(uniprot_ids)
    )

    results <- lapply(uniprot_ids, af_prediction)

    if (anyNA(results)) {
        not_found <- is.na(results)
        n_not_found <- sum(not_found)
        not_found_accessions <- uniprot_ids[not_found]
        if (n_not_found > 6L)
            not_found_accessions <- c(head(not_found_accessions, 5L), "...")
        spdl::info(
            "{} of {} uniprot accessions not found\n  {}",
            n_not_found, length(results),
            paste(sQuote(not_found_accessions, FALSE), collapse = ", ")
        )
        results <- results[!not_found]
    }

    bind_rows(results)
}

af_prediction_view_validate <-
    function(tbl, bfc)
{
    ## AlphaMissenseR 'Suggests: r3dmol'; r3dmol 'Imports: bio3d'
    if (!requireNamespace("r3dmol", quietly = TRUE))
        stop("'af_view_predictions()' requires package 'r3dmol'")

    ## validation
    required_columns <- c(
        "uniprot_id", "protein_variant", "am_pathogenicity", "am_class"
    )
    stopifnot(
        inherits(tbl, "tbl"),
        all(required_columns %in% colnames(tbl)),
        inherits(bfc, "BiocFileCache")
    )
    uniprot_id <-
        distinct(tbl, .data$uniprot_id) |>
        pull("uniprot_id")
    stopifnot(
        `'tbl' must contain a single uniprot_id` =
            identical(length(uniprot_id), 1L)
    )

    uniprot_id
}

#' @rdname AlphaFold
#'
#' @description `af_prediction_view()` summarizes effects of possible
#'     amino acid changes in a single UniProt protein. The changes are
#'     displayed on the AlphaFold-predicted structure.
#'
#' @param tbl
#'
#' A tibble containing information on the UniProt protein and
#' AlphaMissense predicted amino acid effects.
#'
#' For `av_prediction_view()` the tibble must have columns
#' `uniprot_id`, `protein_variant`, `am_pathogenicity`, and
#' `am_class`, as in tibbles returned by `am_data("hg38")` or
#' `am_data("aa_substitutions")`, for instance. The `uniprot_id` must
#' contain a single unique value.
#'
#' For `af_colorfunc_by_position()` the tibble must have columns `pos`
#' and `value`, as described below.
#'
#' @param bfc An object created with `BiocFileCache::BiocFileCache()`,
#'     representing the location used to cache PDB files retrieved by
#'     `av_prediction_view()`. The default is the BiocFileCache
#'     installation-wide location.
#'
#' @details
#'
#' `af_prediction_view()` uses `tbl` to calculate median pathogenicity
#' at each amino acid position, using
#' `am_aa_pathogenicity()`. Predicted protein structure is retrieved
#' from the unique `uniprot_id` using `af_predictions()` and the
#' `pdbUrl` returned by that function. Protein structure is visualized
#' using the r3dmol <https://cran.R-project.org/package=r3dmol>
#' package. Amino acids are colored using `aa_pathogenicity_median`
#' and `af_colorfunc_by_position()` with default palette defined on
#' the interval 0, 1.
#'
#' @return
#'
#' `af_prediction_view()` displays an interactive view of the protein
#' in an RStudio panel or browser tab.
#'
#' @examples
#'
#' ## af_prediction_view()
#'
#' P35557 <-
#'     am_data("hg38") |>
#'     dplyr::filter(uniprot_id == "P35557")
#' af_prediction_view(P35557)
#'
#' ## no AlphaFold prediction for this protein
#' P35555 <-
#'     am_data("aa_substitutions") |>
#'     dplyr::filter(uniprot_id == "P35555")
#' tryCatch({
#'     af_prediction_view(P35555)
#' }, error = identity)
#'
#' @export
af_prediction_view <-
    function(tbl, bfc = BiocFileCache())
{
    uniprot_id <- af_prediction_view_validate(tbl, bfc)

    ## summarize amino acid position pathogenicity
    pathogenicity <- am_aa_pathogenicity(tbl)

    ## retrieve AlphaFold prediction
    prediction <- af_prediction(uniprot_id)
    if (identical(prediction, NA)) {
        stop(spdl::fmt(
            "'af_prediction_view()' could not find UniProt accession '{}'",
            uniprot_id
        ))
    }

    ## retrieve PDB model, using BiocFileCache
    pdb_url <- pull(prediction, "pdbUrl")
    pdb_file <- bfcrpath(bfc, rnames = basename(pdb_url), fpath = pdb_url)
    pdb <- bio3d::read.pdb(pdb_file)

    ## create javascript function for coloring by amino acid position
    colorfunc <-
        pathogenicity |>
        af_colorfunc_by_position(
            "aa_pos", "aa_pathogenicity_median",
            length(pdb$seqres),
            palette_min = 0, palette_max = 1
        )

    ## visualize using r3dmol
    r3dmol::r3dmol() |>
        r3dmol::m_add_model(r3dmol::m_bio3d(pdb)) |>
        r3dmol::m_set_style(
            style = r3dmol::m_style_cartoon(
                arrows = TRUE, colorfunc = colorfunc
            )
        ) |>
        r3dmol::m_zoom_to()
}

#' @rdname AlphaFold
#'
#' @description `af_colorfunc_by_position()` generates a Javascript
#'     function to be used in `rd3mol::m_set_style()` to color
#'     residues by position, e.g., when visualizing median predicted
#'     pathogenicity.
#'
#' @usage
#' af_colorfunc_by_position(
#'     tbl,
#'     pos,
#'     value,
#'     pos_max = NULL,
#'     palette = colorspace::diverging_hcl(11),
#'     palette_min = NULL,
#'     palette_max = NULL
#' )
#'
#' @param pos the symbol or name of the column in `tbl` containing
#'     amino acid residue positions in the protein.
#'
#' @param value the symbol or name of the column in `tbl` containing
#'     values to be used for coloring amino acid residues in the
#'     protein.
#'
#' @param pos_max integer(1) the maximum residue position in the
#'     protein to be visualized. Default: the maximum value in `pos`.
#'
#' @param palette character() vector of colors to be used in
#'     visualization. The default (`colorspace::diverging_hcl(11)`)
#'     produces colors ranging from blue (low) to red (high).
#'
#' @param palette_min numeric(1) the value bounding the minimum
#'     palette color. The default is the minimum of `value`; a common
#'     value when plotting pathogenicity might be `0`.
#'
#' @param palette_max numeric(1) the value bounding the maximum
#'     palette color. The default is the maximum of `value`; a common
#'     value when plotting pathogenicity might be `1`.
#'
#' @details
#'
#' `af_colorfunc_by_position()` uses a template mechanism to inject
#' a vector of position-specific colors into a Javascript function
#' used by `r3dmol::m_set_style()` / `r3dmol::m_style_cartoon()` to
#' color residues by position. Positions for which no color is
#' specified are colored `'gray'`. The template can be seen with
#' `AlphaMissenseR:::js_template("colorfunc")`.
#'
#' @return
#'
#' `af_colorfunc_by_position()` returns a character(1) vector
#' representation of the Javascript function, with color vector
#' injected.
#'
#' @examples
#'
#' ## af_colorfunc_by_position()
#'
#' df <- tibble(
#'     pos = 1 + 1:10, # no color information for position 1
#'     value = 10:1 / 10
#' )
#' colorfunc <- af_colorfunc_by_position(
#'     df,
#'     "pos", "value",
#'     pos_max = 12    # no color information for position 12
#' )
#' cat(colorfunc)
#'
#' ## template used for Javascript function
#' cat(
#'     AlphaMissenseR:::js_template("colorfunc", colors = "..."),
#'     "\n"
#' )
#'
#' @export
af_colorfunc_by_position <-
    function(
        tbl,
        pos, value, pos_max = NULL,
        palette = colorspace::diverging_hcl(11),
        palette_min = NULL,
        palette_max = NULL
    )
{
    ## validation
    stopifnot(
        inherits(tbl, "data.frame"),
        is_scalar_character(pos),
        is_scalar_character(value),
        all(c(pos, value) %in% colnames(tbl)),
        is.null(pos_max) || is_scalar_numeric(pos_max),
        is_character(palette),
        is.null(palette_min) || is_scalar_numeric(palette_min),
        is.null(palette_max) || is_scalar_numeric(palette_max)
    )

    pos <- pull(tbl, pos)
    value <- pull(tbl, value)

    ## default values
    if  (is.null(pos_max))
        pos_max <- max(pos)
    if (is.null(palette_min))
        palette_min <- min(value)
    if (is.null(palette_max))
        palette_max <- max(value)

    ## generate vector of colors
    colors <- character(pos_max)
    colors[] <- "gray" # default
    breaks <- seq(
        palette_min,
        palette_max * (1 + .Machine$double.eps),
        length.out = length(palette) + 1L
    )
    colors[pos] <- palette[findInterval(value, breaks)]

    ## create the function using the template
    colors <- paste(sQuote(colors, q = FALSE), collapse = ", ")
    js_template("colorfunc", colors = colors)
}
