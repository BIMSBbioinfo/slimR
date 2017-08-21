#' searchSLiMs
#'
#' Using regular expressions of SLiMs, find all matches of given SLiM patterns
#' in a given sequence
#'
#' @param sequence Protein sequence using amino acid alphabet
#' @param motifRegex A list object where names of the list are ELM identifiers
#'   and each list item has one ELM regex
#' @param from integer value (default: 1) where to start the search for patterns
#' @param to integer value (default: length of input sequence) where to end the
#'   search
#' @examples
#' data("glutFasta")
#' motifRegex <- list('TRG_ENDOCYTIC_2' = 'Y..[LMVIF]')
#' searchSLiMs(paste(glutFasta), motifRegex)
#' @return A vector of motif hits with the syntax: <MotifIdentifier>:<start
#'   position in sequence>:<end position in sequence> e.g.:
#'   TRG_ENDOCYTIC_2:333:340
#' @export
searchSLiMs <- function(sequence, motifRegex, from = 1, to = nchar(sequence)) {

  seqLen <- nchar(sequence)

  if(from > seqLen | to < 1) {
    stop("Searching coordinates must be within the limit of the sequence.
         Choose values between 1 and ",nchar(sequence),
         " for the 'from' and 'to' arguments respectively\n")
  }

  if(from < 1) {
    warning("Starting position to look for regex matches must be >= 1 ... setting the value of 'from' to 1")
    from <- 1
  }

  if(to > seqLen) {
    warning("Ending position to look for regex matches must be <= length of the input sequence.
            Setting the value of 'to' to ",nchar(sequence),"\n")
    to <- seqLen
  }

  sequence <- substr(sequence, from, to)
  if (from > 1) {
    #to avoid matching N terminal motifs when the searched
    #string is not a prefix of the original sequence
    sequence <- paste0("XXX", sequence)
    from <- from - 3
  }

  if(to < seqLen) {
    #to avoid matching C terminal motifs when the searched
    #string is not a suffix of the original sequence
    sequence <- paste0(sequence, "XXX")
  }

  hits <- do.call(rbind, lapply(X = names(motifRegex), function(x) {
     m <- slimR::locateAllRegex(sequence = sequence,
                         pattern = motifRegex[[x]])
     if(!is.null(m)) {
       m$start <- m$start + from - 1
       m$end <- m$end + from - 1
       m$SLiM <- x
       return(m)
      }
     }))
  return(hits)
}

#' locateAllRegex
#'
#' Finds the positions and substrings of all overlapping
#' matches of regexes in a given sequence
#'
#' @param sequence A character vector
#' @param pattern PERL like regular expression
#' @return data.table data.frame object
#' @importFrom data.table data.table
#' @export
locateAllRegex <- function (sequence, pattern) {
  startPos <- pracma::refindall(s = sequence, pat = pattern)
  matchSeq <- unlist(lapply(startPos, function(x) {
    pracma::regexp(s = substring(sequence, x), pat = pattern, once = TRUE)$match
  }))
  endPos <- startPos + nchar(matchSeq) - 1
  if(!is.null(startPos)) {
    return(data.table::data.table('start' = startPos, 'end' = endPos, 'match' = matchSeq))
  }
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
#' @param variants A data.frame consisting of minimum four columns:
#'   1.uniprotAccession, 2.wtAA, 3.mutAA, 4.pos where pos is the mutation position in
#'   the sequence, wtAA is the wild-type amino acid (one letter code) in the
#'   sequence and mutAA is the mutant amino acid (one letter code).
#' @return data.table data.frame object
#' @importFrom dplyr setdiff
#' @examples
#' c <- findMotifChanges(sequence = paste(glutFasta),
#'                       variants = glutMutations,
#'                      motifRegex = motifRegex)
#' @export
findMotifChanges <- function(sequence, variants, motifRegex = slimR::motifRegex) {

  if(sum(c("uniprotAccession", "wtAA", "mutAA", "pos") %in% colnames(variants)) != 4) {
    stop("Seems like the variants data.frame does not contain all of the required columns:
         uniprotAccession, wtAA, mutAA, and pos")
  }

  variants$uniprotAccession <- as.character(variants$uniprotAccession)
  variants$wtAA <- as.character(variants$wtAA)
  variants$mutAA <- as.character(variants$mutAA)
  variants$pos <- as.numeric(variants$pos)

  wtMotifsAll <- slimR::searchSLiMs(sequence = sequence, motifRegex = motifRegex)
  if(!is.null(wtMotifsAll)) {
    wtMotifsAll$ID <- apply(wtMotifsAll, 1,
                            function(x) gsub(' ', '', paste0(c(x['SLiM'], x['start'], x['end']), collapse = ':')))
  }

  #for each variant find the list of motif hits that would be gained/lost
  #if the sequence had that substitution mutation
  change <- do.call(rbind, lapply(1:nrow(variants), function(i) {
    #cat("i:",i,"\n")
    wtAA <- variants$wtAA[i]
    mutAA <- variants$mutAA[i]
    pos <- variants$pos[i]

    mutSeq <- slimR::mutateSequence(sequence = sequence,
                                    pos = pos,
                                    wtAA = wtAA,
                                    mutAA = mutAA)

    mutMotifs <- slimR::searchSLiMs(sequence = mutSeq, motifRegex = motifRegex, from = pos - 51, to = pos + 50)
    #filter mutMotifs for those that overlap variant position

    if(!is.null(mutMotifs)) {
      mutMotifs <- mutMotifs[pos >= start & pos <= end]
      mutMotifs$ID <- apply(mutMotifs, 1,
                            function(x) gsub(' ', '', paste0(c(x['SLiM'], x['start'], x['end']), collapse = ':')))
    }

    #filter wtMotifs for those that overlap variant position
    if(!is.null(wtMotifsAll)) {
      wtMotifs <- wtMotifsAll[pos >= start & pos <= end]
    }

    gained <- data.table()
    lost <- data.table()

    if(is.null(mutMotifs) & is.null(wtMotifsAll)) {
      return(NULL)
    } else if (is.null(mutMotifs)) {
      lost <- wtMotifs
    } else if (is.null(wtMotifs)) {
      gained <- mutMotifs
    } else {
      lost <- wtMotifs[which(wtMotifs$ID %in% setdiff(wtMotifs$ID, mutMotifs$ID)),]
      gained <- mutMotifs[which(mutMotifs$ID %in% setdiff(mutMotifs$ID, wtMotifs$ID)),]
    }

    if(nrow(lost) == 0 & nrow(gained) == 0){
      result <- data.table()
    } else if (nrow(lost) == 0) {
      gained$type <- 'gained'
      result <- gained
    } else if (nrow(gained) == 0) {
      lost$type <- 'lost'
      result <- lost
    } else {
      lost$type <- 'lost'
      gained$type <- 'gained'
      result <- rbind(lost,gained)
    }
    if(nrow(result) > 0){
      result$wtAA <- wtAA
      result$mutAA <- mutAA
      result$pos <- pos
      result$uniprotAccession <- variants$uniprotAccession[i]
      return(result)
    }
  }))

  if(!is.null(change)){
    change <- subset(change, select = c('uniprotAccession', 'wtAA', 'mutAA', 'pos',
                                        'SLiM', 'start', 'end', 'type', 'match'))

    colnames(change) <- c('uniprotAccession', 'wtAA', 'mutAA', 'pos',
                          'SLiM', 'SLiM_start', 'SLiM_end',
                          'change', 'SLiM_Sequence')

    change$RegEx <- paste(motifRegex[change$SLiM])

    unChanged <- dplyr::setdiff(x = variants[,c('uniprotAccession', 'wtAA', 'mutAA', 'pos')],
                                y = change[,c('uniprotAccession', 'wtAA', 'mutAA', 'pos')])
    if(nrow(unChanged) > 0) {
      unChanged <- cbind(unChanged, 'SLiM' = 'None', 'SLiM_start' = 0, 'SLiM_end' = 0,
                         'change' = 'NoChange', 'SLiM_Sequence' = 'None', 'RegEx' = 'None')
      return(rbind(change, unChanged))
    } else {
      return(change)
    }
  } else {
    unChanged <- cbind(variants[,c('uniprotAccession', 'wtAA', 'mutAA', 'pos')],
                       'SLiM' = 'None', 'SLiM_start' = 0, 'SLiM_end' = 0,
                       'change' = 'NoChange', 'SLiM_Sequence' = 'None', 'RegEx' = 'None')
    return(unChanged)
  }
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
