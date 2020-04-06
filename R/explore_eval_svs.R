##' Opens a Shiny app with a dynamic table that contains FP, FN and TP SVs. Clicking on a
##' SV (row in the table) generates a simplified representation of the variants in the region.
##' Useful to explore the calls, the truth and what was called TP/FN/FP.
##' @title Interactive exploration of FP/FN/TP SVs
##' @param eval.o output from svevalOl.R
##' @param ucsc.genome the genome version for the UCSC Genome Browser automated link.
##' @return Starts a Shiny app in a web browser.
##' @author Jean Monlong
##' @importFrom magrittr %>%
##' @importFrom rlang .data
##' @export
explore_eval_svs <- function(eval.o, ucsc.genome='hg38'){

  svs.gr = lapply(eval.o$svs, function(esvtype){
    lapply(names(esvtype), function(etype){
      gr = esvtype[[etype]]
      gr$eval = etype
      gr
    })
  })
  names(svs.gr) = NULL
  svs.gr = do.call(what=c, args=unlist(svs.gr))
  
  ## List for plot_ranges
  svs.grl = GenomicRanges::split(svs.gr, svs.gr$eval)

  ## data.frame for the dynamic table
  svs.dfl = lapply(svs.grl, function(svs.gr){
    svs.gr %>% as.data.frame %>%
      dplyr::mutate(chr=.data$seqnames,
                    coord=paste0(.data$chr,':',.data$start,'-',.data$end)) %>%
      dplyr::arrange(.data$chr, .data$start) %>% 
      dplyr::select(.data$coord, .data$type, .data$size) %>%
      dplyr::mutate(type=factor(.data$type))
  })
  names(svs.dfl) = names(svs.grl)
    
  ui <- shiny::fluidPage(
                 shiny::titlePanel("Eval output exploration"),
                 shiny::sidebarLayout(
                          shiny::sidebarPanel(
                                   width=3,
                                   shiny::numericInput('flank', 'Flanks (bp):', 500, min=1, max=1e4, step=50),
                                   shiny::radioButtons('evalt', 'Eval type:', c('FP','FN','TP','TP.baseline'), 'FP'),
                                   DT::dataTableOutput('vartable')),
                          shiny::mainPanel(
                                   width=9,
                                   shiny::htmlOutput('url', class='btn btn-default action-button shiny-bound-input'),
                                   shiny::fluidRow(shiny::plotOutput('ranges'))
                                 )
                        )
               )
                 
  server <- function(input, output) {
    output$vartable <- DT::renderDataTable(
                             svs.dfl[[input$evalt]],
                             filter='top',
                             rownames=FALSE,
                             options=list(pageLength=15),
                             selection='single')
    output$ranges <- shiny::renderPlot({
      sv.sel = svs.dfl[[input$evalt]][input$vartable_rows_selected,]
      ggp = NULL
      if(!is.na(sv.sel$coord[1])){
        gr = GenomicRanges::GRanges(sv.sel$coord[1])
        ggp = plot_ranges(svs.grl, gr, maxgap=input$flank, pt.size=3, lab.size=6) +
          ggplot2::theme(text=ggplot2::element_text(size=18)) +
          ggplot2::scale_colour_brewer(palette='Set1')
      }
      return(ggp)
    })
    output$url = shiny::renderText({
      sv.sel = svs.dfl[[input$evalt]][input$vartable_rows_selected,]
      coord = sv.sel$coord[1]
      if(!is.na(coord)){
        gr = GenomicRanges::GRanges(coord)
        chr.name = as.character(GenomicRanges::seqnames(gr))[1]
        ggp = plot_ranges(svs.grl, gr, maxgap=input$flank, pt.size=3, lab.size=6) +
          ggplot2::theme(text=ggplot2::element_text(size=18))
        x.lims = round(ggp$scales$scales[[1]]$limits)
        coord = paste0(chr.name, ':', max(0, x.lims[1]), '-', x.lims[2])
      }
      paste0('<a href="https://genome.ucsc.edu/cgi-bin/hgTracks?db=', ucsc.genome,
             '&position=', coord, '" target="_blank">UCSC Genome Browser</a>')
    })
  }

  shiny::shinyApp(ui=ui, server=server)
}
