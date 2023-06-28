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

  ## extract informatin and make GRanges from the eval.o input list
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

  ## data.frame for the dynamic table (for each type of variant: FP, FN, TP, TP.baseline)
  svs.dfl = lapply(svs.grl, function(svs.gr){
    svs.gr %>% as.data.frame %>%
      dplyr::mutate(chr=.data$seqnames,
                    coord=paste0(.data$chr,':',.data$start,'-',.data$end)) %>%
      dplyr::arrange(.data$chr, .data$start) %>% 
      dplyr::select(.data$coord, .data$type, .data$size) %>%
      dplyr::mutate(type=factor(.data$type)) %>%
      dplyr::sample_frac(1)
  })
  names(svs.dfl) = names(svs.grl)

  ## color palette
  col.v = c(FP="#e41a1c", FN="#984ea3", TP="#377eb8", TP.baseline="#4daf4a")
  ## shapes for each SV type
  shape.v = c(DEL=15, INS=17, DUP=18, INV=16, BND=4)
  if(!all(names(eval.o$svs) %in% names(shape.v))){
    shape.v = seq(1, length(eval.o$svs))
    names(shape.v) = names(eval.o$svs)
  }
  
  ## app interface
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

  ## server side of the app
  server <- function(input, output) {
    ## output dynamic table
    output$vartable <- DT::renderDataTable(
                             svs.dfl[[input$evalt]],
                             filter='top',
                             rownames=FALSE,
                             options=list(pageLength=15),
                             selection='single')
    ## ggplot representing the SV of interest in the region
    output$ranges <- shiny::renderPlot({
      sv.sel = svs.dfl[[input$evalt]][input$vartable_rows_selected,]
      ggp = NULL
      if(!is.na(sv.sel$coord[1])){
        gr = GenomicRanges::GRanges(sv.sel$coord[1])
        ggp = plot_ranges(svs.grl, gr, maxgap=input$flank, pt.size=3, lab.size=6) +
          ggplot2::theme(text=ggplot2::element_text(size=18)) +
          ggplot2::scale_colour_manual(values=col.v) + 
          ggplot2::scale_shape_manual(values=shape.v)
      }
      return(ggp)
    })
    ## link to the region in UCSC browser. becomes a button above using class definition
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

  ## launch app
  shiny::shinyApp(ui=ui, server=server)
}
