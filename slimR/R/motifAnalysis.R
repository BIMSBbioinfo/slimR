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
searchSLiMs <- function(sequence, motifRegex) {
  hits <- lapply(X = motifRegex, function(x) {locateAllRegex(sequence = sequence,
                                                             pattern = x)})
  hits <- hits[lapply(hits, nrow)  > 0]
  unlist(lapply(X = c(1:length(names(hits))),
                FUN = function(x) {paste0(names(hits)[x],
                                          ':', hits[[x]]$start,
                                          ':', hits[[x]]$end)}))
}

#' locateAllRegex
#'
#' Finds the positions and substrings of all overlapping
#' matches of regexes in a given sequence
#'
#' @param sequence A character vector
#' @param pattern PERL like regular expression
#' @export
locateAllRegex <- function (sequence, pattern) {
  startPos <- pracma::refindall(s = sequence, pat = pattern)
  matchSeq <- unlist(lapply(startPos, function(x) {
    pracma::regexp(s = substring(sequence, x), pat = pattern, once = TRUE)$match
    }))
  endPos <- startPos + nchar(matchSeq) - 1
  return(data.frame('start' = startPos, 'end' = endPos, 'match' = matchSeq, stringsAsFactors = FALSE))
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


#' findMotifChanges
#'
#' Find out which SLiMs are gained or lost (no longer matching the regex
#' pattern) via point amino acid substitutions in protein sequences
#'
#' @param sequence A character string of amino acid sequence
#' @param variants A data.frame consisting of minimum three columns: 1.wtAA,
#'   2.mutAA, 3.pos where pos is the mutation position in the sequence, wtAA is
#'   the wild-type amino acid (one letter code) in the sequence and mutAA is the
#'   mutant amino acid (one letter code).
#' @examples
#' c <- findMotifChanges(sequence = paste(glutFasta),
#'                       variants = glutMutations,
#'                      motifRegex = motifRegex)
#' @export
findMotifChanges <- function(sequence, variants, motifRegex) {

  #find SLiMs in the wild-type sequence
  wtMotifs <- searchSLiMs(sequence = sequence, motifRegex = motifRegex)
  #for each variant find the list of motif hits that would be gained/lost
  #if the sequence had that substitution mutation

  lostMotifs <- c()
  gainedMotifs <- c()

  for (i in 1:nrow(variants)) {
    wtAA <- as.character(variants$wtAA[i])
    mutAA <- as.character(variants$mutAA[i])
    pos <- as.integer(variants$pos[i])

    mutSeq <- mutateSequence(sequence = sequence,
                             pos = pos,
                             wtAA = wtAA,
                             mutAA = mutAA)

    mutMotifs <- searchSLiMs(sequence = mutSeq, motifRegex = motifRegex)

    lost <- setdiff(wtMotifs, mutMotifs)
    if (length(lost) > 0) {
      lostMotifs <- c(lostMotifs,
                      paste(wtAA, mutAA, pos,
                            lost, 'lost', sep = ':'))
    }

    gained <- setdiff(mutMotifs, wtMotifs)

    if (length(gained) > 0) {
      gainedMotifs <- c(gainedMotifs,
                        paste(wtAA, mutAA, pos,
                              gained, 'gained', sep = ':'))
    }
  }
  change <- c(lostMotifs, gainedMotifs)
  if (length(change) > 0) {
    wtSeq <- unlist(strsplit(sequence, split = ''))
    change <- data.frame(do.call(rbind,
                                 strsplit(change, ':')),
                         stringsAsFactors = FALSE)
    colnames(change) <- c('wtAA', 'mutAA', 'pos',
                          'SLiM', 'SLiM_start', 'SLiM_end',
                          'change')
    change$RegEx <- motifRegex[change$SLiM]
    change$SLiM_Sequence <- lapply(X = 1:nrow(change),
                                   FUN = function (x) {
                                     paste(wtSeq[change$SLiM_start[x]:change$SLiM_end[x]],
                                           collapse = '')})
    return(change)
  }
}

#' findMotifChangesBulk
#'
#' Find out which SLiMs are gained or lost (no longer matching the regex
#' pattern) via point amino acid substitutions in protein sequences
#'
#' @param variants A Granges object of variants parsed from Humsavar
#' using getHumSavar() function or a subset of it with the same structure
#'   consisting of minimum two meta-data columns: 1.wtAA,
#'   2.mutAA, wtAA is the wild-type amino acid (one letter code) in the
#'   sequence and mutAA is the mutant amino acid (one letter code).
#' @param uniprotDataDir The folder that is used to download/store uniprot
#' data files
#' @export
findMotifChangesBulk <- function (variants, uniprotDataDir) {

  ptm <- proc.time()

  motifChanges <- list()

  #load regular expressions for motifs
  data("motifRegex")

  uniprotAccessions <- unique(as.vector(seqnames(variants)))
  downloadUniprotFiles(uniprotAccessions = uniprotAccessions, outDir = uniprotDataDir, format = 'fasta')

  pb <- txtProgressBar(min = 0, max = length(uniprotAccessions), style = 3)

  for (i in 1:length(uniprotAccessions)) {
    setTxtProgressBar(pb, i)
    uniAcc <- uniprotAccessions[i]
    fastaFile <- file.path(uniprotDataDir, 'fasta', paste0(uniAcc, '.fasta'))

    if(file.exists(fastaFile)) {
      sequence <- paste(Biostrings::readAAStringSet(filepath = fastaFile, format = 'fasta'))

      myVariants <- variants[seqnames(variants) == uniAcc]

      df <- data.frame('wtAA' = as.vector(myVariants$wtAA),
                       'mutAA' = as.vector(myVariants$mutAA),
                       'pos' = start(myVariants))

      changes <- findMotifChanges(sequence = sequence, variants = df, motifRegex = motifRegex)
      if (!is.null(changes)) {
        changes$uniprotAcc <- uniAcc
        motifChanges[[length(motifChanges)+1]] <- changes
        names(motifChanges)[length(motifChanges)] <- uniAcc
      }
    }
  }

  close(pb)
  proc.time() - ptm
  return(motifChanges)
}

