##' Opens a Shiny app with a dynamic table that contains input SVs. Clicking on a
##' SV (row in the table) generates a simplified representation of the variation
##' graph around this SV. The number of flanking nodes (context) can be increased
##' if necessary, e.g. for large insertions. vg needs to be installed
##' (https://github.com/vgteam/vg).
##' @title Interactive exploration of SVs in a variation graph
##' @param svs either a GRanges with SVs (e.g. from \code{readSVvcf}) or the path to
##' a VCF file.
##' @param xg the path to the xg object of the variation graph.
##' @param ucsc.genome the genome version for the UCSC Genome Browser automated link.
##' @return Starts a Shniy app in a web browser.
##' @author Jean Monlong
##' @importFrom magrittr %>%
##' @importFrom rlang .data
##' @export
ivg_sv <- function(svs, xg, ucsc.genome='hg38'){

  if(is.character(svs) & length(svs)==1){
    message("Reading VCF file...")
    svs = readSVvcf(svs)
  }
  
  svs.df = as.data.frame(svs) %>%
    dplyr::mutate(chr=.data$seqnames,
                  coord=paste0(.data$chr,':',.data$start,'-',.data$end)) %>% 
    dplyr::filter(.data$chr=='chr10') %>% 
    dplyr::select(.data$coord, .data$type, .data$size)

  ui <- shiny::fluidPage(
    shiny::titlePanel("IVG-SV"),
    shiny::sidebarLayout(
      shiny::sidebarPanel(width=4,
        shiny::numericInput('context', 'Context (nodes):', 3, min=1, max=100, step=1),
        DT::dataTableOutput('vartable')),
      shiny::mainPanel(width=8,
               shiny::column(3, shiny::selectInput('zoom', 'Zoom:',
                                                   c('100%','200%', '400%'), '100%')),
        shiny::column(3, shiny::htmlOutput('url', class='btn btn-default action-button shiny-bound-input')),
        shiny::uiOutput('svg')
      )
    )
  )

  server <- function(input, output) {
    output$vartable <- DT::renderDataTable(
                             svs.df,
                             filter='top',
                             rownames=FALSE,
                             options=list(pageLength=15),
                             selection='single')
    output$diagram <- DiagrammeR::renderGrViz({
      sv.sel = svs.df[input$vartable_rows_selected,]
      find.o = system2('vg', args=c('find', '-x', xg, '-p', sv.sel$coord[1],
                                    '-c', input$context), stdout='temp.vg')
      mod.o = system2('vg', args=c('mod', '-u', 'temp.vg'), stdout='temp2.vg')
      view.o = system2('vg', args=c('view', '-Sdp', 'temp2.vg'), stdout=TRUE)
      DiagrammeR::grViz(view.o)
    })
    output$svg = shiny::renderUI({
      DiagrammeR::grVizOutput('diagram', width=input$zoom)
    })
    output$url = shiny::renderText({
      sv.sel = svs.df[input$vartable_rows_selected,]      
      paste0('<a href="https://genome.ucsc.edu/cgi-bin/hgTracks?db=', ucsc.genome,
             '&position=', sv.sel$coord[1], '" target="_blank">UCSC Genome Browser</a>')
    })
  }

  shiny::shinyApp(ui=ui, server=server)
}
