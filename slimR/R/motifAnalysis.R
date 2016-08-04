#' searchMotifs
#'
#' Using regular expressions of SLiMs, find all matches of given SLiM patterns
#' in a given sequence
#'
#' @param sequence Protein sequence using amino acid alphabet
#' @param motifRegex A list object where names of the list are ELM identifiers
#'  and each list item has one ELM regex
#' @examples
#' data(sampleFasta)
#' motifRegex <- list('TRG_ENDOCYTIC_2' = 'Y..[LMVIF]')
#' #or
#' #motifRegex <- as.list(as.vector(elms$RegEx))
#' #names(motifRegex) <- as.vector(elms$ELM_Identifier)
#' searchSLiMs(sampleFasta, motifRegex)
#' @return A list of data.frame objects.
#' @export
searchSLiMs <- function (sequence, motifRegex) {
  hits <- lapply(X = motifRegex,
                 FUN = function(x) { stringr::str_locate_all(sequence, x)})
  hits <- lapply(X = hits, FUN = function (x) { df <- data.frame(x[[1]]) })
  hits <- hits[lapply(hits, nrow)  > 0]
  #for( i in 1:length(hits)) {
   # hits[[i]]$motif <- names(hits)[i]
  #}
  #hitsDF <- do.call(rbind, hits)
  #rownames(hitsDF) <- c(1:nrow(hitsDF))

  #return(hitsDF)
  return(hits)
}



