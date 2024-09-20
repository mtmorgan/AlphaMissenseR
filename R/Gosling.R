gosling_view_bar <-
    function(rng, track_data)
{
    ## define single track
    track_bar <- shiny.gosling::add_single_track(
        width = 800,
        height = 180,
        data = track_data,
        mark = "bar",
        x = shiny.gosling::visual_channel_x(
            field = "start", type = "genomic", axis = "bottom"
        ),
        xe = shiny.gosling::visual_channel_x(field = "end", type = "genomic"),
        y = shiny.gosling::visual_channel_y(
            field = "am_pathogenicity", type = "quantitative", axis = "right"
        ),
        color = shiny.gosling::visual_channel_color(
            field = "am_pathogenicity", type = "quantitative"
        ),
        tooltip = shiny.gosling::visual_channel_tooltips(
            shiny.gosling::visual_channel_tooltip(
                field = "REF", type = "nominal", alt = "Reference"
            ),
            shiny.gosling::visual_channel_tooltip(
                field = "ALT", type = "nominal", alt = "Alternative / Mutation"
            ),
            shiny.gosling::visual_channel_tooltip(
                field = "am_pathogenicity", type = "quantitative",
                alt = "AM_Pathogenicity Score", format = "0.2"
            )
        ),
        size = list(value = 5)
    )

    shiny.gosling::compose_view(
        layout = "linear",
        xDomain = list(
            chromosome = as.character(GenomicRanges::seqnames(rng)),
            interval = c(GenomicRanges::start(rng), GenomicRanges::end(rng))
        ),
        tracks = track_bar
    )
}

gosling_view_lollipop <-
    function(rng, track_data)
{
    ## Define categories and color mapping
    categories <- c("likely_benign", "ambiguous", "likely_pathogenic")
    colormapping <- c("#89d5f5", "gray", "#f56c6c")

    ## Define multi tracks
    track_ref <- shiny.gosling::add_single_track(
        data = track_data,
        mark = "rect",
        x = shiny.gosling::visual_channel_x(
            field = "start", type = "genomic", axis = "top"
        ),
        xe = shiny.gosling::visual_channel_x(
            field = "end", type = "genomic"
        ),
        size = list(value = 50),
        stroke = "lightgrey",
        strokeWidth = list(value = 1),
        opacity = list(value = 0.3)
    )

    track_alt <- shiny.gosling::add_single_track(
        data = track_data,
        mark = "point",
        x = shiny.gosling::visual_channel_x(
            field = "start", type = "genomic", axis = "top"
        ),
        xe = shiny.gosling::visual_channel_x(field = "end", type = "genomic"),
        y = shiny.gosling::visual_channel_y(
            field = "am_class", type="nominal", domain= categories,
            axis = "left",baseline = "ambiguous"
        ),
        text = list(field = "ALT", type = "nominal"),
        size = list(value = 5),
        tooltip = shiny.gosling::visual_channel_tooltips(
            shiny.gosling::visual_channel_tooltip(
                field = "REF", type = "nominal", alt = "Reference"
            ),
            shiny.gosling::visual_channel_tooltip(
                field = "ALT", type = "nominal", alt = "Alternative / Mutation"
            ),
            shiny.gosling::visual_channel_tooltip(
                field = "am_pathogenicity",
                type = "quantitative",
                alt = "AM_Pathogenicity Score",
                format = "0.2"
            )
        )
    )

    ## Compose view
    shiny.gosling::compose_view(
        width = 800,
        height = 180,
        multi = TRUE,
        layout = "linear",
        xDomain = list(
            chromosome = as.character(GenomicRanges::seqnames(rng)),
            interval = c(GenomicRanges::start(rng), GenomicRanges::end(rng))
        ),
        alignment = "overlay",
        color = shiny.gosling::visual_channel_color(
            field = "am_class",
            type = "nominal",
            domain = categories,
            baseline = "ambiguous",
            range = colormapping,
            legend = TRUE
        ),
        tracks = shiny.gosling::add_multi_tracks(track_ref, track_alt)
    )
}

gosling_app <-
    function(view)
{
    ## Create Shiny app
    ui <- shiny::fluidPage(
        shiny.gosling::use_gosling(clear_files = FALSE),
        shiny.gosling::goslingOutput("gosling_plot")
    )

    server <- function(input, output, session) {
        output$gosling_plot <- shiny.gosling::renderGosling({
            shiny.gosling::gosling(component_id = "component_3", view)
        })
    }

    ## Return the Shiny app
    shiny::shinyApp(ui, server)
}

#' @rdname Gosling
#'
#' @title Plot a GPos or GRanges Object Using Shiny and Gosling
#'
#' @description `plot_granges()` creates a Shiny app that displays a
#'     Gosling plot of a given GRanges object.  It visualizes genomic
#'     ranges with both rectangle and point representations, and
#'     allows for customization of the plot title and subtitle.
#'
#' @details
#'
#' The function supports 2 types of plots as selected through the
#' `plot_type` argument: (1) `"bar"`: a barplot-like view of the
#' pathogenicity score at each position, similar to a sequencing
#' coverage plot, and (2) `"lollipop"`: a lollipop plot focusing on
#' the pathogenicity classification (ambiguous, benign, pathogenic) at
#' each position. It requires that the [`GenomicRanges::GRanges`]
#' object has the following metadata columns: 'am_class' (effect
#' classification),'am_pathogenicity' (pathogenicity score), 'ALT'
#' (alternative allele) and 'REF' (reference allele).
#'
#' @param granges A [`GenomicRanges::GRanges`] (e.g.,
#'     [`GenomicRanges::GPos`]) object containing the genomic ranges
#'     to be plotted.
#'
#' @param title character(1) Title of the plot. Default is
#'     "GRanges Plot".
#'
#' @param subtitle character(1) The subtitle of the plot. Default is
#'     "Stacked nucleotide example".
#'
#' @param plot_type character(1) Select the type of gosling plot.
#'     Default is "bar"; available types are described in Details.
#'
#' @param create_app logical(1) Produce a Shiny app for Gosling
#'     visualization. Default `TRUE` when used interactively.
#'
#' @return `gosling_plot()` with `create_app = TRUE` returns a
#'     `shinyApp` object that, when run, displays the Gosling
#'     plot. When `create_app = FALSE`, the (invisible) return value
#'     represents the gosling object to be used in Shiny apps.
#'
#' @note This function requires the `shiny`, `shiny.gosling`, and
#'     `GenomicRanges` packages to be installed.
#'
#' @examples
#' ## Create a sample GRanges object from AlphamissenseR
#' gpos <-
#'     am_data("hg38") |>
#'     filter(uniprot_id == "Q1W6H9") |>
#'     to_GPos()
#'
#' ## Plot `gpos` using shiny.gosling
#' if (requireNamespace("shiny.gosling", quietly = TRUE)) {
#'    gosling_plot(
#'        gpos, title = "Q1W6H_track", subtitle = "bar plot example",
#'        plot_type = "bar", create_app = FALSE
#'    )
#' }
#'
#' @importFrom tools R_user_dir
#'
#' @importFrom methods as
#'
#' @export
gosling_plot <-
    function(
        granges,
        title = "GRanges Plot",
        subtitle = "Stacked nucleotide example",
        plot_type = c("bar", "lollipop"),
        create_app = interactive()
    )
{
    ## Validate input
    stopifnot(
        is(granges, "GRanges"),
        isScalarCharacter(title),
        isScalarCharacter(subtitle),
        isScalarLogical(create_app),
        requireNamespace("GenomicRanges", quietly = TRUE),
        requireNamespace("shiny.gosling", quietly = TRUE)
    )
    plot_type <- match.arg(plot_type)

    ## Turns out gr must be coerced to GRanges(), will look into this later
    granges <- as(granges, "GRanges")
    ## Get range from GRanges object
    rng <- range(granges)

    ## avoid writing '.gosling' to current working directory; awkward
    ## because shiny.gosling uses the working directory.
    cache <- R_user_dir("AlphaMissenseR", which = "cache")
    gosling_cache <- file.path(cache, ".gosling")
    if (!dir.exists(gosling_cache)) {
        status <- dir.create(gosling_cache, recursive = TRUE)
        if (!status)
            stop("failed to create gosling cache '", gosling_cache, "'")
    }
    original_wd <- getwd()
    on.exit(setwd(original_wd))
    setwd(cache) # so that '.gosling' is available to shiny.gosling

    ## Prepare track data
    track_data <- shiny.gosling::track_data_gr(
        granges,
        chromosomeField = "seqnames",
        genomicFields = c("start", "end")
    )

    composed_view <- switch(
        plot_type,
        bar = gosling_view_bar(rng, track_data),
        lollipop = gosling_view_lollipop(rng, track_data),
        stop("unknown `plot_type = '", plot_type, "'`")
    )

    ## Arrange into view
    arranged_view <- shiny.gosling::arrange_views(
        title = title,
        subtitle = subtitle,
        views = composed_view
    )

    if (create_app) {
        gosling_app(arranged_view)
    } else {
        ## return arranged_view invisibly
        invisible(arranged_view)
    }
}
