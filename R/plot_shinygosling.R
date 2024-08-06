#' @rdname plot_shinygosling
#'
#' @title Plot a GRanges object using Shiny and Gosling
#'
#' @description This function creates a Shiny app that displays a
#'     Gosling plot of a given GRanges object.  It visualizes genomic
#'     ranges with both rectangle and point representations, and
#'     allows for customization of the plot title and subtitle.

#'     The function creates two tracks: a reference track with
#'     rectangles and an alternative track with points. It assumes
#'     that the GRanges object has 'am_class' and 'ALT' metadata
#'     columns. The 'am_class' column is used for coloring the points,
#'     while 'ALT' is used for the text of the points.

#' @param gr A GRanges object containing the genomic ranges to
#'     be plotted.
#' @param title Character string. The title of the plot. Default is
#'     "GRanges Plot".
#' @param subtitle Character string. The subtitle of the plot. Default
#'     is "Stacked nucleotide example".
#'
#' @return A Shiny app object that, when run, displays the Gosling
#'     plot.
#'
#' @note This function requires the shiny, shiny.gosling, and GenomicRanges
#'     packages to be installed.
#'
#' @examples
#' if (requireNamespace("GenomicRanges")) {
#'
#' ## Create a sample GRanges object from AlphamissenseR
#' gpos <-
#'    am_data("hg38") |>
#'    filter(uniprot_id == "Q1W6H9") |>
#'    to_GPos()
#'
#' ## Plot the GRanges object
#' plot_granges(gpos, mode = "bar", title = "Q1W6H9 track", subtitle 
#' = "bar plot example")
#' }
#'
#'
#' @export
plot_granges <-
    function(gr,
             title = "GRanges Plot",
             subtitle = "Stacked nucleotide example")
{
    ## Validate input
    stopifnot(
        "Input must be a GRanges or GPos object" = 
            inherits(gr, c("GRanges", "GPos"))
    )
    stopifnot(
        "Option must be a " = 
            class(title)=="character"
        )
    stopifnot(
        "Option must be a " = 
            class(subtitle)=="character"
    )
    
    
            
    
    ## Define categories and color mapping
    categories <- c("likely_benign", "ambiguous", "likely_pathogenic")
    colormapping <- c("#89d5f5", "gray", "#f56c6c")

    # Get range from GRanges object
    r <- range(gr)

    # This fixes the bug if .gosling directory does not already exist
    if (!dir.exists(".gosling")){
        dir.create(".gosling")
    }
    
    # Prepare track data
    track_data <- track_data_gr(
        gr,
        chromosomeField = "seqnames",
        genomicFields = c("start", "end")
    )

    # Define tracks
    track_ref <- add_single_track(
        data = track_data,
        mark = "rect",
        x = visual_channel_x(field = "start", type = "genomic", axis = "bottom"),
        xe = visual_channel_x(field = "end", type = "genomic"),
        size = list(value = 50),
        stroke = "lightgrey",
        strokeWidth = list(value = 1),
        opacity = list(value = 0.3)
    )

    track_alt <- add_single_track(
        data = track_data,
        mark = "point",
        x = visual_channel_x(field = "start", type = "genomic", axis = "bottom"),
        xe = visual_channel_x(field = "end", type = "genomic"),
        y = visual_channel_y(field = "am_class", type = "nominal", axis = "right"),
        text = list(field = "ALT", type = "nominal"),
        size = list(value = 5)
    )

    # Compose view
    composed_view_a <- compose_view(
        width = 800,
        height = 180,
        multi = TRUE,
        layout = "linear",
        xDomain = list(
            chromosome = as.character(seqnames(r)),
            interval = c(start(r), end(r))
        ),
        alignment = "overlay",
        color = visual_channel_color(
            field = "am_class",
            type = "nominal",
            domain = categories,
            range = colormapping,
            legend = TRUE
        ),
        tracks = add_multi_tracks(track_ref, track_alt)
    )

    # Arrange view
    arranged_view3 <- arrange_views(
        title = title,
        subtitle = subtitle,
        views = composed_view_a
    )

    # Create Shiny app
    ui <- fluidPage(
        use_gosling(clear_files = FALSE),
        goslingOutput("gosling_plot")
    )

    server <- function(input, output, session) {
        output$gosling_plot <- renderGosling({
            gosling(
                component_id = "component_3",
                arranged_view3
            )
        })
    }

    # Return the Shiny app
    shinyApp(ui, server)
}
