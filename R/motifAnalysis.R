#' searchSLiMs
#'
#' Using regular expressions of SLiMs, find all matches of given SLiM patterns
#' in a given sequence
#'
#' @param sequence Protein sequence using amino acid alphabet
#' @param motifRegex A list object where names of the list are ELM identifiers
#'   and each list item has one ELM regex
#' @examples
#' data("glutFasta")
#' motifRegex <- list('TRG_ENDOCYTIC_2' = 'Y..[LMVIF]')
#' searchSLiMs(paste(glutFasta), motifRegex)
#' @return A vector of motif hits with the syntax:
#' <MotifIdentifier>:<start position in sequence>:<end position in sequence>
#' e.g.: TRG_ENDOCYTIC_2:333:340
#' @export
searchSLiMs <- function(sequence, motifRegex) {
  hits <- lapply(X = motifRegex, function(x) {locateAllRegex(sequence = sequence,
                                                             pattern = x)})
  hits <- hits[lapply(hits, nrow)  > 0]
  if (length(hits) > 0) {
    unlist(lapply(X = c(1:length(names(hits))),
                  FUN = function(x) {paste0(names(hits)[x],
                                            ':', hits[[x]]$start,
                                            ':', hits[[x]]$end)}))
  } else {
    return(c())
  }
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
findMotifChanges <- function(sequence, variants, motifRegex = slimR::motifRegex) {

  if(sum(c("wtAA", "mutAA", "pos") %in% colnames(variants)) != 3) {
    stop("Seems like the variants data.frame does not contain all of the required columns:
         wtAA, mutAA, and pos")
  }

  #find SLiMs in the wild-type sequence
  wtMotifs <- searchSLiMs(sequence = sequence, motifRegex = motifRegex)
  #for each variant find the list of motif hits that would be gained/lost
  #if the sequence had that substitution mutation

  lostMotifs <- c()
  gainedMotifs <- c()
  neutralVars <- c() #Variants that don't cause any changes in motif content

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

    if(length(lost) == 0 & length(gained) == 0){
      neutralVars <- c(neutralVars, paste(wtAA, mutAA, pos,
                             "None:0:0",
                             'NoChange', sep = ':'))
    }

  }
  change <- c(lostMotifs, gainedMotifs, neutralVars)

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
                                   if(change$change[x] == 'NoChange'){
                                     "None"
                                   } else {
                                     paste(wtSeq[change$SLiM_start[x]:change$SLiM_end[x]],
                                           collapse = '')
                                   }
                                   })
  change$SLiM_start <- as.numeric(change$SLiM_start)
  change$SLiM_end <- as.numeric(change$SLiM_end)
  change$pos <- as.numeric(change$pos)
  return(change)
}

#' findMotifChangesMulti
#'
#' A wrapper function that runs findMotifChanges function for multiple inputs
#' using multiple cores
#'
#' @param sequences List of strings where the strings are amino-acid sequences
#'   and the list names are uniprot accessions
#' @param motifRegex List of slim regular expressions. By default, the built-in
#'   slimR::motifRegex from the ELM database is used.
#' @param variants A data.frame consisting of minimum four columns: 1.
#'   uniprotAccession 2.wtAA, 3.mutAA, 4.pos where pos is the mutation position
#'   in the sequence, wtAA is the wild-type amino acid (one letter code) in the
#'   sequence and mutAA is the mutant amino acid (one letter code).
#' @param nodeN Number of cores needed to run the analysis (default: 1)
#' @return List of data.frame objects. One data.frame per each uniprot accession
#' @examples
#' sequences <- downloadUniprotFiles(uniprotAccessions = c('P04637', 'P11166'),
#'                          format = 'fasta', overwrite = 'FALSE', nodeN = 2)
#' variants <- slimR::glutMutations
#' motifRegex <- list('mymotif' = '.ABSDBASDBAS.')
#' motifChanges <- findMotifChangesMulti(sequences = sequences,
#'                                      variants = variants,
#'                                      motifRegex = motifRegex,
#'                                      nodeN = 1)
#' @export
findMotifChangesMulti <- function(sequences,
                                  variants,
                                  motifRegex = slimR::motifRegex,
                                  nodeN = 1) {

  if(sum(c("uniprotAccession", "wtAA", "mutAA", "pos") %in% colnames(variants)) != 4) {
    stop("Seems like the variants data.frame does not contain all of the required columns:
         uniprotAccession, wtAA, mutAA, and pos")
  }

  cl <- parallel::makeCluster(nodeN)
  doParallel::registerDoParallel(cl)

  motifChanges <- foreach (i=1:length(sequences), .inorder = TRUE) %dopar% {
    uni <- names(sequences)[i]
    if (uni %in% variants$uniprotAccession) {
      result <- findMotifChanges(sequence = sequences[[i]],
                                  variants = variants[variants$uniprotAccession == uni,],
                                  motifRegex = motifRegex
                                  )
    } else {
      result <- 'No variants found'
    }
    result
  }
  names(motifChanges) <- names(sequences)

  stopCluster(cl)
  return(motifChanges)
}
