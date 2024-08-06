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
#' @param plot_type Character string. Select the type of gosling plot. 
#' Default is "bar"
#'     - "bars": Stacked bar plot with height based on pathogenicity score
#'     - "lolipop": variation of a bar chart where the bar is replaced with a 
#'     line and a dot at the end to show mutation variations.
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
             subtitle = "Stacked nucleotide example",
             plot_type = "bars")
{
    ## Validate input
    stopifnot(
        "Input must be a GRanges or GPos object" = 
            inherits(gr, c("GRanges", "GPos"))
    )
    stopifnot(
        "Title must be a scalar string" = 
            is.character(title) && length (title) == 1 && !is.na(title) && 
            nzchar(title)
        )
    stopifnot(
        "Subtitle must be a scalar string" = 
            is.character(subtitle) && length (subtitle) == 1 && !is.na(subtitle)
            && nzchar(subtitle)
    )
    stopifnot(
        "Type must be a scalar string" = 
            is.character(plot_type) && length (plot_type) == 1 && !is.na(plot_type)
        && nzchar(plot_type)
    )
    
    ## Define categories and color mapping
    categories <- c("likely_benign", "ambiguous", "likely_pathogenic")
    colormapping <- c("#89d5f5", "gray", "#f56c6c")

    ## Get range from GRanges object
    r <- range(gr)
    
    ## Prepare track data
    track_data <- track_data_gr(
        gr,
        chromosomeField = "seqnames",
        genomicFields = c("start", "end")
    )

    ## This fixes the bug if .gosling directory does not already exist
    cache_dir <- file.path(tools::R_user_dir("AlphaMissenseR", which = "cache"), ".gosling")
    if (!dir.exists(cache_dir))
        ## TODO: check return value to ensure directory is created successfully
        dir.create(cache_dir, recursive = TRUE)
    
    ## trigger the option for bars or lolipop        
    if (plot_type =="bars"){
        #define single track
        track_bar <- add_single_track(
            width = 800,
            height = 180,
            data = track_data,
            mark = "bar",
            x = visual_channel_x(
                field = "start", type = "genomic", axis = "bottom"
            ),
            xe = visual_channel_x(field = "end", type = "genomic"),
            y = visual_channel_y(
                field = "am_pathogenicity", type = "quantitative", axis = "right"
            ),
            color = visual_channel_color(
                field = "am_pathogenicity",
                type = "quantitative"
            ),
            tooltip = visual_channel_tooltips(
                visual_channel_tooltip(field = "REF", type = "nominal",
                                       alt = "Reference"),
                visual_channel_tooltip(field = "ALT", type = "nominal",
                                       alt = "Alternative / Mutation"),
                visual_channel_tooltip(
                    field = "am_pathogenicity",
                    type = "quantitative",
                    alt = "AM_Pathogenicity Score",
                    format = "0.2"
                ) ),
            size = list(value = 5)
        )
        
        composed_view <- compose_view(
            layout = "linear",
            xDomain = list(chromosome = as.character(seqnames(r)),
                           interval = c(start(r), end(r))),
            tracks = track_bar
            
        )
        
    ## other track    
    } else if (plot_type == "lolipop"){

    ## Define multi tracks
    track_ref <- add_single_track(
        data = track_data,
        mark = "rect",
        x = visual_channel_x(field = "start", type = "genomic", axis = "top"),
        xe = visual_channel_x(field = "end", type = "genomic"),
        size = list(value = 50),
        stroke = "lightgrey",
        strokeWidth = list(value = 1),
        opacity = list(value = 0.3)
    )

    track_alt <- add_single_track(
        data = track_data,
        mark = "point",
        x = visual_channel_x(field = "start", type = "genomic", axis = "top"),
        xe = visual_channel_x(field = "end", type = "genomic"),
        y = visual_channel_y(field = "am_class", type="nominal", 
                       domain= categories, axis = "left",baseline = "ambiguous" ),
        text = list(field = "ALT", type = "nominal"),
        size = list(value = 5),
        tooltip = visual_channel_tooltips(
            visual_channel_tooltip(field = "REF", type = "nominal",
                                   alt = "Reference"),
            visual_channel_tooltip(field = "ALT", type = "nominal",
                                   alt = "Alternative / Mutation"),
            visual_channel_tooltip(
                field = "am_pathogenicity",
                type = "quantitative",
                alt = "AM_Pathogenicity Score",
                format = "0.2"
            ) )
    )

    ## Compose view
    composed_view <- compose_view(
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
            baseline = "ambiguous",
            range = colormapping,
            legend = TRUE
        ),
        tracks = add_multi_tracks(track_ref, track_alt)
    )
    
    } 
    ## trigger check
    else {
        stop("Invalid plot_type. Use 'bars' or 'lolipop'")
    }
    ## Arrange into view
    arranged_view3 <- arrange_views(
        title = title,
        subtitle = subtitle,
        views = composed_view
    )

    ## Create Shiny app
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

    ## Return the Shiny app
    shinyApp(ui, server)
}
