##' Opens a Shiny app with a dynamic table that contains SV clusters. Clicking on a
##' SV (row in the table) generates a simplified representation of the SV clusters variants in the region.
##' Useful to explore the clusters and understand if under-merging.
##' @title Interactive exploration of FP/FN/TP SVs
##' @param svs SVs clustered by site with clusterSVs
##' @param annot optional: a GRanges object with regions to highlight (e.g. repeats)
##' @param ucsc.genome the genome version for the UCSC Genome Browser automated link.
##' @param graph.height height of the graph in the shiny app, in pixel. Default is 400
##' @param ucsc.url the URL/mirror for the UCSC Genome Browser
##' @param subtitle A subtitle to add. Default is NULL
##' @return Starts a Shiny app in a web browser.
##' @author Jean Monlong
##' @importFrom magrittr %>%
##' @importFrom rlang .data
##' @export
explore_sv_clusters <- function(svs, annot=NULL, ucsc.genome='hg38', graph.height=400,
                                ucsc.url='genome.ucsc.edu', subtitle=NULL){

  svs.df = svs %>% as.data.frame %>% 
  dplyr::group_by(.data$seqnames, .data$svsite, .data$type, .data$clique) %>%
    dplyr::summarize(ac=length(.data$size),
                     start=min(.data$start),
                     end=max(.data$end),
                     size.min=min(.data$size),
                     size.max=max(.data$size),
                     size=round(mean(.data$size)),
                     .groups='drop')
  svs.gr = GenomicRanges::makeGRangesFromDataFrame(svs.df, keep.extra.columns=TRUE)
  
  ## data.frame for the dynamic table (for each type of variant: FP, FN, TP, TP.baseline)
  svs.df = svs.df %>% 
      dplyr::mutate(coord=paste0(.data$seqnames,':',.data$start,'-',.data$end)) %>%
      dplyr::arrange(.data$seqnames, .data$start) %>% 
      dplyr::select(.data$coord, .data$type, .data$size, .data$ac, .data$clique) %>%
      dplyr::mutate(type=factor(.data$type), clique=factor(.data$clique)) %>%
      dplyr::sample_frac(1)
  
  ## shapes for each SV type
  shape.v = c(DEL=15, INS=17, DUP=18, INV=16, BND=4)
  if(!all(unique(svs.df$type) %in% names(shape.v))){
    shape.v = c(15, 17, 18, 16, 4, seq(1, length(unique(svs.df$type))))[1:length(unique(svs.df$type))]
    names(shape.v) = unique(svs.df$type)
  }

  svtab.cols = c('seqnames', 'start', 'end', 'svsite', 'type',
                 'size', 'ac', 'clique')

  p.title = "SV clustering exploration"
  if(!is.null(subtitle)){
    p.title = paste("SV clustering exploration - ", subtitle)
  }
  
  ## app interface
  ui <- shiny::fluidPage(
                 shiny::titlePanel(p.title),
                 shiny::sidebarLayout(
                          shiny::sidebarPanel(
                                   width=3,
                                   shiny::numericInput('flank', 'Flanks (bp):', 500, min=1, max=1e4, step=50),
                                   shiny::numericInput('minsize', 'Minimum SV size (bp):', 30, min=10, max=Inf, step=1),
                                   shiny::numericInput('minac', 'Minimum number of allels:', 0, min=0, max=Inf, step=1),
                                   DT::dataTableOutput('vartable')),
                          shiny::mainPanel(
                                   width=9,
                                   shiny::htmlOutput('url', class='btn btn-default action-button shiny-bound-input'),
                                   shiny::fluidRow(shiny::plotOutput('ranges', height=graph.height)),
                                   DT::dataTableOutput('svtable')
                                 )
                        )
               )

  ## server side of the app
  server <- function(input, output) {
    getSVdf = reactive({
      svs.df[which(svs.df$size >= input$minsize & svs.df$ac >= input$minac),]
    })

    ## output dynamic table
    output$vartable <- DT::renderDataTable(
                             getSVdf(),
                             filter='top',
                             rownames=FALSE,
                             options=list(pageLength=15),
                             selection='single')

    getSVgr = reactive({
      sv.sel = getSVdf()[input$vartable_rows_selected,]
      gr.svs = GRanges()
      if(!is.na(sv.sel$coord[1])){
        gr = GenomicRanges::GRanges(sv.sel$coord[1])
        gr.svs = svs.gr[which(svs.gr$size >= input$minsize & svs.gr$ac >= input$minac)]
        gr.svs = IRanges::subsetByOverlaps(gr.svs, gr, maxgap=input$flank)
      }
      gr.svs
    })

    getSVgr.fortable = reactive({
      gr.svs = getSVgr()
      df = data.frame()
      if(length(gr.svs) > 0){
        df = subset(svs, svsite %in% gr.svs$svsite) %>%
          as.data.frame %>%
          dplyr::mutate(type=factor(.data$type),
                        clique=factor(.data$clique),
                        svsite=factor(.data$svsite)) %>%
          arrange(.data$svsite, .data$size)
        df = df[, svtab.cols]
      }
      df
    })
    
    ## ggplot representing the SV of interest in the region
    output$ranges <- shiny::renderPlot({
      sv.sel = getSVdf()[input$vartable_rows_selected,]
      gr.svs = getSVgr()
      ggp = NULL
      if(length(gr.svs)>0){
        gr = GenomicRanges::GRanges(sv.sel$coord[1])
        grl = GenomicRanges::split(gr.svs, gr.svs$svsite)
        ggp = plot_ranges(grl, gr, maxgap=input$flank, pt.size=3, lab.size=6,
                          gr.annot=annot) +
          ggplot2::theme(text=ggplot2::element_text(size=18)) +
          ggplot2::scale_shape_manual(values=shape.v)
      }
      return(ggp)
    })
    ## link to the region in UCSC browser. becomes a button above using class definition
    output$url = shiny::renderText({
      sv.sel = getSVdf()[input$vartable_rows_selected,]
      gr.svs = getSVgr()
      coord = ''
      if(length(gr.svs)>0){
        gr = GenomicRanges::GRanges(sv.sel$coord[1])
        grl = GenomicRanges::split(gr.svs, gr.svs$svsite)
        chr.name = as.character(GenomicRanges::seqnames(gr))[1]
        ggp = plot_ranges(grl, gr, maxgap=input$flank, pt.size=3, lab.size=6,
                          gr.annot=annot) +
          ggplot2::theme(text=ggplot2::element_text(size=18))
        x.lims = round(ggp$scales$scales[[1]]$limits)
        coord = paste0(chr.name, ':', max(0, x.lims[1]), '-', x.lims[2])
      }
      paste0('<a href="https://', ucsc.url, '/cgi-bin/hgTracks?db=', ucsc.genome,
             '&position=', coord, '" target="_blank">UCSC Genome Browser</a>')
    })

    ## output dynamic table
    output$svtable <- DT::renderDataTable(
                            getSVgr.fortable(),
                            filter='top',
                            rownames=FALSE,
                            options=list(pageLength=15),
                            selection='single')

  }

  ## launch app
  shiny::shinyApp(ui=ui, server=server)
}
