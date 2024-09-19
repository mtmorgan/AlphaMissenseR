#' @rdname plot_granges
#'
#' @title Plot a GRanges object using Shiny and Gosling
#'
#' @description This function creates a Shiny app that displays a
#'     Gosling plot of a given GRanges object.  It visualizes genomic
#'     ranges with both rectangle and point representations, and
#'     allows for customization of the plot title and subtitle.

#' @details
#' The function supports 2 types of plots as selected through the `plot_type` 
#' argument: (1) a barplot-like view of the pathogenicity score at each position, 
#' similar to a sequencing coverage plot, and (2) a lollipop plot focusing on 
#' the pathogenicity classification (ambiguous, benign, pathogenic) at each 
#' position. It requires that the [`GRanges`] object has the following metadata
#' columns: 'am_class' (effect classification),'am_pathogenicity' 
#' (pathogenicity score), 'ALT' (alternative allele) and 'REF' 
#' (reference allele).

#' @param gr A [`GRanges`] object containing the genomic ranges to
#'     be plotted.
#' @param title character(1) Title of the plot. Default is
#'     "GRanges Plot".
#' @param subtitle character(1) The subtitle of the plot. Default
#'     is "Stacked nucleotide example".
#' @param plot_type character(1) Select the type of gosling plot. 
#' Default is "bars".
#'     - "bars": Stacked bar plot with height based on pathogenicity score,
#'     - "lollipop": variation of a bar chart where the bar is replaced with a 
#'     line and a dot at the end to show mutation variations.
#'
#' @return A `shinyApp` object that, when run, displays the Gosling
#'     plot.
#'
#' @note This function requires the `shiny`, `shiny.gosling`, and `GenomicRanges`
#'     packages to be installed.
#'
#' @examples
#' if (requireNamespace("GenomicRanges")) {
#'
#'    ## Create a sample GRanges object from AlphamissenseR
#'    gpos <-
#'        am_data("hg38") |>
#'        filter(uniprot_id == "Q1W6H9") |>
#'        to_GPos()
#'
#'    ## Plot the GRanges object
#'     plot_granges(
#'        gpos, mode = "bar", title = "Q1W6H_track",
#'        subtitle = "bar plot example"
#'     )
#'}
#'
#'
#'
#' @export
plot_granges <-
    function(gr_input,
             title = "GRanges Plot",
             subtitle = "Stacked nucleotide example",
             plot_type = "bars")
{
    ## Validate input
    stopifnot(
            is(granges, "GRanges"),
            isScalarCharacter(title),
            isScalarCharacter(subtitle),
            isScalarCharacter(plot_type)
        )
    
    ## Define categories and color mapping
    categories <- c("likely_benign", "ambiguous", "likely_pathogenic")
    colormapping <- c("#89d5f5", "gray", "#f56c6c")
    
    ## Turns out gr must be coerced to GRanges(), will look into this later
    gr_input <- as(gr_input, "GRanges")
    ## Get range from GRanges object
    r <- range(gr_input)
    
    # This fixes the bug if .gosling directory does not already exist
    if (!dir.exists(".gosling")){
        dir.create(".gosling")
    }

    ## Prepare track data
    track_data <- shiny.gosling::track_data_gr(
        gr_input,
        chromosomeField = "seqnames",
        genomicFields = c("start", "end")
    )

    
    ## trigger the option for bars or lollipop        
    if (identical(plot_type, "bars")){
        #define single track
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
                field = "am_pathogenicity",
                type = "quantitative"
            ),
            tooltip = shiny.gosling::visual_channel_tooltips(
                shiny.gosling::visual_channel_tooltip(field = "REF", type = "nominal",
                                       alt = "Reference"),
                shiny.gosling::visual_channel_tooltip(field = "ALT", type = "nominal",
                                       alt = "Alternative / Mutation"),
                shiny.gosling::visual_channel_tooltip(
                    field = "am_pathogenicity",
                    type = "quantitative",
                    alt = "AM_Pathogenicity Score",
                    format = "0.2"
                ) ),
            size = list(value = 5)
        )
        
        composed_view <- shiny.gosling::compose_view(
            layout = "linear",
            xDomain = list(chromosome = as.character(seqnames(r)),
                           interval = c(start(r), end(r))),
            tracks = track_bar
            
        )
        
    ## other option    
    } else if (identical(plot_type, "lollipop")){

    ## Define multi tracks
    track_ref <- shiny.gosling::add_single_track(
        data = track_data,
        mark = "rect",
        x = shiny.gosling::visual_channel_x(field = "start", type = "genomic", axis = "top"),
        xe = shiny.gosling::visual_channel_x(field = "end", type = "genomic"),
        size = list(value = 50),
        stroke = "lightgrey",
        strokeWidth = list(value = 1),
        opacity = list(value = 0.3)
    )

    track_alt <- shiny.gosling::add_single_track(
        data = track_data,
        mark = "point",
        x = shiny.gosling::visual_channel_x(field = "start", type = "genomic", axis = "top"),
        xe = shiny.gosling::visual_channel_x(field = "end", type = "genomic"),
        y = shiny.gosling::visual_channel_y(field = "am_class", type="nominal", 
                       domain= categories, axis = "left",baseline = "ambiguous" ),
        text = list(field = "ALT", type = "nominal"),
        size = list(value = 5),
        tooltip = shiny.gosling::visual_channel_tooltips(
            shiny.gosling::visual_channel_tooltip(field = "REF", type = "nominal",
                                   alt = "Reference"),
            shiny.gosling::visual_channel_tooltip(field = "ALT", type = "nominal",
                                   alt = "Alternative / Mutation"),
            shiny.gosling::visual_channel_tooltip(
                field = "am_pathogenicity",
                type = "quantitative",
                alt = "AM_Pathogenicity Score",
                format = "0.2"
            ) )
    )

    ## Compose view
    composed_view <- shiny.gosling::compose_view(
        width = 800,
        height = 180,
        multi = TRUE,
        layout = "linear",
        xDomain = list(
            chromosome = as.character(seqnames(r)),
            interval = c(start(r), end(r))
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
    ## trigger check
    else {
        stop("Invalid plot_type. Use 'bars' or 'lollipop'")
    }
    ## Arrange into view
    arranged_view3 <- shiny.gosling::arrange_views(
        title = title,
        subtitle = subtitle,
        views = composed_view
    )

    ## Create Shiny app
    ui <- shiny::fluidPage(
        shiny.gosling::use_gosling(clear_files = FALSE),
        shiny.gosling::goslingOutput("gosling_plot")
    )

    server <- function(input, output, session) {
        output$gosling_plot <- shiny.gosling::renderGosling({
            shiny.gosling::gosling(
                component_id = "component_3",
                arranged_view3
            )
        })
    }

    ## Return the Shiny app
    shiny::shinyApp(ui, server)
}
