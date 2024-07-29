#' @rdname plot_shinygosling
#' Plot a GRanges object using Shiny and Gosling
#'
#' This function creates a Shiny app that displays a Gosling plot of a given GRanges object.
#' It visualizes genomic ranges with both rectangle and point representations, and allows
#' for customization of the plot title and subtitle.
#'
#' @param as_GRanges A GRanges object containing the genomic ranges to be plotted.
#' @param title Character string. The title of the plot. Default is "GRanges Plot".
#' @param subtitle Character string. The subtitle of the plot. Default is "Stacked nucleotide example".
#'
#' @return A Shiny app object that, when run, displays the Gosling plot.
#'
#' @details The function creates two tracks: a reference track with rectangles and an
#' alternative track with points. It assumes that the GRanges object has 'am_class' and 'ALT'
#' metadata columns. The 'am_class' column is used for coloring the points, while 'ALT' is
#' used for the text of the points.
#'
#' @note This function requires the shiny, gosling, and GenomicRanges packages to be installed.
#'
#' @examples
#' \dontrun{
#' library(GenomicRanges)
#' 
#' # Create a sample GRanges object
#' gr <- GRanges(
#'   seqnames = c("chr1", "chr2", "chr1"),
#'   ranges = IRanges(start = c(1, 3, 5), end = c(2, 4, 6)),
#'   am_class = c("likely_benign", "ambiguous", "likely_pathogenic"),
#'   ALT = c("A", "T", "G")
#' )
#' 
#' # Plot the GRanges object
#' plot_granges(gr, title = "My GRanges Plot", subtitle = "Custom subtitle")
#' }
#'
#' @import shiny
#' @import gosling
#' @import GenomicRanges
#'
#' @export
plot_granges <- function(as_GRanges, title = "GRanges Plot", subtitle = "Stacked nucleotide example") {
    # Define categories and color mapping
    categories <- c("likely_benign", "ambiguous", "likely_pathogenic")
    colormapping <- c("#029F73", "gray", "#CB3B8C")
    
    # Get range from GRanges object
    get_range <- range(as_GRanges)
    
    # Prepare track data
    track_data <- track_data_gr(
        as_GRanges,
        chromosomeField = "seqnames",
        genomicFields = c("start", "end")
    )
    
    # Define tracks
    track_ref <- add_single_track(
        data = track3_data,
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
            chromosome = as.character(get_range@seqnames@values),
            interval = c(get_range@ranges@start, get_range@ranges@start + get_range@ranges@width)
        ),
        alignment = "overlay",
        color = visual_channel_color(
            field = "am_class",
            type = "nominal",
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