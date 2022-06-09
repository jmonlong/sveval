##' Cluster SVs based on overlap/similarity
##'
##' SVs are overlapped with each other. A graph is then built where nodes (SVs) are connected is
##' they overlap/match. A cluster is a component in this graph.
##'
##' To reduce the memory usage, the SVs are first grossly clustered into batches. The actual clustering
##' (and graph construction) is then performed separately for each batch, potentially in parallel.
##'
##' @title Cluster SVs based on overlap/similarity
##' @param svs.gr A GRanges with SVs (for example read by \code{readSVvcf} or \code{readSVvcf.multisamps})
##' @param min.rol minimum reciprocal overlap for deletions and other "ranges" SVs. Default is 0.1
##' @param max.ins.dist maximum distance for insertions to be clustered.
##' @param range.seq.comp compare sequence instead of overlapping deletions/inversion/etc. Default is FALSE.
##' @param ins.seq.comp compare sequence instead of insertion sizes. Default is FALSE.
##' @param nb.cores number of processors to use. Default is 1.
##' @param simprep optional simple repeat annotation. Default is NULL. If non-NULL, GRanges to be used to
##' @param batch.maxsize batch size to aim. To reduce memory usage, see Details. Default is 500.
##' @param log.level the level of information in the log. Default is "CRITICAL" (basically no log).
##' @return the svs.gr object annotated with two new columns:
##' \item{svsite}{the ID of the cluster (or SV site)}
##' \item{clique}{is this cluster a clique, i.e. all SVs overlapping/matching all other SVs in the cluster}
##' @author Jean Monlong
##' @export
##' @importFrom magrittr %>%
##' @importFrom rlang .data
##' @examples
##' \dontrun{
##'
##' svs = readSVvcf('svs.vcf.gz', keep.ids=TRUE)
##' svs = clusterSVs(svs)
##' 
##' }
clusterSVs <- function(svs.gr, min.rol=.8, max.ins.dist=20,
                       range.seq.comp=FALSE, ins.seq.comp=FALSE, simprep=NULL,
                       nb.cores=1, batch.maxsize = 500,
                       log.level=c('CRITICAL', 'WARNING', 'INFO')){
  logging::setLevel(log.level[1])

  ## are all the SVs connected/overlapping/matched to every other SVs in the cluster?
  is_clique <- function(amat){
    diag(amat) = 1
    all(amat==1)
  }

  ## Function to overlap and annotate a batch of SVs
  ## Used on slices of the full dataset to reduce memory consumption
  annotateSVs <- function(gr){
    ## overlaps SVs in the set
    ol.df = svOverlap(gr, gr, min.ol=min.rol, min.del.rol=min.rol,
                      max.ins.dist=max.ins.dist, 
                      range.seq.comp=range.seq.comp, ins.seq.comp=ins.seq.comp,
                      simprep=simprep,
                      log.level=log.level) %>%
      as.data.frame %>%
      dplyr::filter(.data$queryHits!=.data$subjectHits)
    ## make graph
    adj.mat = Matrix::sparseMatrix(ol.df$queryHits, ol.df$subjectHits, x=rep(1, nrow(ol.df)), dims=c(length(gr), length(gr)))
    rownames(adj.mat) = colnames(adj.mat) = gr$svid
    ## extract components
    cmp.o = igraph::components(igraph::graph_from_adjacency_matrix(adj.mat, mode='undirected'))
    ## annotate GRanges, checking if components are cliques
    svid.idx = 1:length(gr)
    names(svid.idx) = as.character(gr$svid)
    gr.cl = lapply(unique(as.numeric(cmp.o$membership)), function(cmp){
      cmp.ii = which(cmp.o$membership == cmp)
      gr.ss = gr[svid.idx[rownames(adj.mat)[cmp.ii]]]
      gr.ss$svsite = rownames(adj.mat)[cmp.ii[1]]
      if(length(cmp.ii) < 3){
        gr.ss$clique = TRUE
      } else {
        gr.ss$clique = is_clique(adj.mat[cmp.ii, cmp.ii])
      }
      return(gr.ss)
    })
    return(do.call(c, gr.cl))
  }

  ## add SV ids is missing
  if(is.null(svs.gr$svid)){
    svs.gr$svid = paste0('sv_', 1:length(svs.gr))
  }
  if(any(duplicated(svs.gr$svid))){
    stop('Input SVs not suitable for clustering: duplicated SV ids')
  }
  if(any(is.na(svs.gr$svid))){
    stop('Input SVs not suitable for clustering: some SV ids are missing')
  }
  
  ## quick clustering to split dataset for analysis
  ## make batches of the SV clusters to have {batch.maxsize} variants
  ## (possibly many more if a cluster can't be split)
  cl.gr = GenomicRanges::reduce(svs.gr, min.gapwidth=2*max.ins.dist)
  cl.gr$n = GenomicRanges::countOverlaps(cl.gr, svs.gr)
  cl.gr = cl.gr[order(cl.gr$n)]
  cl.gr$batch = cumsum(cl.gr$n) %>% cut(seq(0,sum(cl.gr$n)+batch.maxsize, batch.maxsize))
  svs.per.batch = unlist(tapply(cl.gr$n, cl.gr$batch, sum))
  logging::loginfo(paste('SV clustering. Number of batches: ', length(svs.per.batch)))
  logging::loginfo(paste('SV clustering. Average number of SVs per batch: ', mean(svs.per.batch)))
  logging::loginfo(paste('SV clustering. Maximum number of SVs per batch: ', max(svs.per.batch)))

  ## merge each batch of large SV cluster
  svs.m = parallel::mclapply(unique(cl.gr$batch), function(bb){
    cl.gr = cl.gr[which(cl.gr$batch==bb)]
    annotateSVs(IRanges::subsetByOverlaps(svs.gr, cl.gr))
  }, mc.cores=nb.cores)

  ## return the bound results
  return(do.call(c, svs.m))
}
