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

  unlist(lapply(X = c(1:length(names(hits))),
                          FUN = function(x) {paste0(names(hits)[x],
                                                    ':', hits[[x]]$start,
                                                    ':', hits[[x]]$end)}))
}

#' mutateSequence
#'
#' Given a protein sequence, create a mutated copy of the sequence at the
#' desired position and amino acids.
#' @param sequence Character string of amino acids
#' @param pos The amino acid number (positive integer) at which a mutation
#'   should be created
#' @param wtAA One letter code of the Wild-type amino-acid (to double check if
#'   the desired mutation matches the actual coordinates of the wild-type
#'   sequence)
#' @param mutAA One letter code of the mutant amino acid
#' @return A character string where wtAA is replaced with mutAA if the
#'   amino-acid at the given position of the given sequence actually has the
#'   same amino acid.
#' @export
mutateSequence <- function (sequence, pos, wtAA, mutAA) {
  if ( pos < 1 | pos > nchar(sequence)) {
    warning ('Mutation position must be between 1 and ',nchar(sequence),'\n')
  } else {
    wtSeq <- unlist(strsplit(sequence, split = ''))
    if (wtSeq[pos] != wtAA) {
      warning ('WT aminoacid at position',pos,
            'does not match the actual amino acid in the given sequence',
            wtSeq[pos],'\n')
    } else {
      mutSeq <- wtSeq
      mutSeq[pos] <- mutAA
      return(paste(mutSeq,collapse = ''))
    }
  }
}


#' findMotifKnockIns
#'
#' Find out which SLiMs are created de novo via point amino acid substitutions
#' in protein sequences
#'
#' @param sequence A character string of amino acid sequence
#' @param variants A data.frame consisting of minimum three columns: 1.wtAA,
#'   2.mutAA, 3.pos where pos is the mutation position in the sequence, wtAA is
#'   the wild-type amino acid in the sequence and mutAA is the mutant amino
#'   acid.
#' @export
findMotifChanges <- function(sequence, variants, motifRegex) {

  #find SLiMs in the wild-type sequence
  wtMotifs <- searchSLiMs(sequence = sequence, motifRegex = motifRegex)
  #for each variant find the list of motif hits that would be gained/lost
  #if the sequence had that substitution mutation
  for (i in 1:length(variants)) {
    wtAA <- as.character(variants$wtAA[i])
    mutAA <- as.character(variants$mutAA[i])
    pos <- as.integer(variants$pos[i])
    mutSeq <- mutateSequence(sequence = sequence,
                             pos = pos,
                             wtAA = wtAA,
                             mutAA = mutAA)
    mutMotifs <- searchSLiMs(sequence = mutSeq, motifRegex = motifRegex)
    lostMotifs <- setdiff(wtMotifs, mutMotifs)
    gainedMotifs <- setdiff(mutMotifs, wtMotifs)

    print(variants[i,])
    cat('Lost Motifs:',lostMotifs,'\n')
    cat('Gained Motifs:',gainedMotifs,'\n')
  }
}

#elms <- getElmClasses()
#motifRegex <- as.list(as.vector(elms$RegEx))
#names(motifRegex) <- as.vector(elms$ELM_Identifier)

#vars <- getHumSavar()
#variants <- vars[vars$uniprotAcc == 'P11166',]

findMotifChanges(sequence = paste(sampleFasta), variants = variants, motifRegex = motifRegex)


