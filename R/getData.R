#' parseMutation
#'
#' Given a vector of mutation substitutions
#' (e.g. "p.His160Arg") -> split "p.His160Arg" into "H 160 R"
#' @param mutations Vector of mutations
#' @return A data.frame object where columns are wtAA, mutAA, and pos
#' @examples
#' mut <- c('p.His160Arg', 'p.Pro12Leu', 'p.Pro14Ter')
#' mut.df <- parseMutation(mut)
#' @export
parseMutation <- function (mutations) {
  mutations <- gsub(pattern = '^p.', replacement = '', x = mutations)
  pos <- stringr::str_match(pattern = '\\d+', string = mutations)
  df <- data.frame(do.call(rbind, stringr::str_split(pattern = '\\d+', string = mutations)))
  df$pos <- as.numeric(pos)
  colnames(df) <- c('wtAA', 'mutAA', 'pos')
  df$wtAA <- slimR::aaTable$oneLetterCode[match(df$wtAA, slimR::aaTable$threeLetterCode, nomatch = )]
  df$mutAA <- slimR::aaTable$oneLetterCode[match(df$mutAA, slimR::aaTable$threeLetterCode)]
  return(df)
}


#' getHumSavar
#'
#' Download and parse protein mutation data from UniProt
#' @param outdir Path to the folder, where the downloaded file should be saved.
#' @return A Granges object containing the coordinates of mutated sites in
#'   proteins
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom stringr str_split
#' @importFrom utils download.file
#' @return A GRanges object where seqnames are Uniprot sequence accessions and
#'   ranges are the coordinates of the variants.
#' @examples
#' \dontrun{
#' hs <- getHumSavar(getwd())
#' }
#' @export
getHumSavar <- function (outdir) {
  variantFile <- file.path(outdir, 'humsavar.txt')
  if (!file.exists(variantFile)) {
    download.file(url = 'www.uniprot.org/docs/humsavar.txt',
                  destfile = variantFile)
  } else {
    warning("humsavar.txt exists at the target location: ",variantFile,
            ", a new one won't be downloaded. Remove the existing
            file and re-run the function to update the file")
  }

  #skip first 50 lines which don't contain mutation data
  dat <- readLines(con = variantFile)[-(1:50)]

  #grep the lines with relevant variant data
  mut <- grep(pattern = "Polymorphism|Disease|Unclassified",
              x = dat,
              perl = TRUE,
              value = TRUE)
  mut <- data.frame(
    do.call(rbind,
            stringr::str_split(string = mut,
                               n = 7,
                               pattern = '\\s+')))

  colnames(mut) <- c('geneName', 'uniprotAccession',
                     'FTId', 'change',
                     'variant', 'dbSNP', 'diseaseName')

  parsedMut <- parseMutation(mutations = mut$change)

  mut <- cbind(mut, parsedMut)

  #some of the mutations may not have been parsed correctly (due to errors in
  #humsavar file). So, some amino acids may have NA values after conversion.
  #such entries are althogether deleted in the merged data frmae)
  mut <- mut[!(is.na(mut$wtAA) | is.na(mut$mutAA) | is.na(mut$pos)),]

  mut <- GenomicRanges::makeGRangesFromDataFrame(df = mut,
                                                 keep.extra.columns = TRUE,
                                                 ignore.strand = TRUE,
                                                 seqnames.field = 'uniprotAccession',
                                                 start.field = 'pos', end.field = 'pos')
  return(mut)
}


### Functions to download and map ClinVar data to Uniprot sequences
#' getClinVarData
#'
#' This function will fetch the clinvar data from the given url and parse the
#' contents of the downloaded file.
#'
#' @param url The url to the ftp location of the clinvar dataset (tab delimited
#'   variant summary file) (e.g.
#'   ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz)
#'
#' @param overwrite default: FALSE, boolean value to decide if a fresh download
#'   should overwrite the existing file
#' @return Path to the downloaded file
#' @importFrom utils download.file
#' @examples
#' \dontrun{
#' #download
#' filePath <- getClinVarData()
#' #read first 100 lines
#' cv_head <- data.table::fread(text = readLines(filePath, n = 100))
#' }
#' @export
getClinVarData <- function(url = paste0('ftp://ftp.ncbi.nlm.nih.gov/',
                                        'pub/clinvar/tab_delimited/',
                                        'variant_summary.txt.gz'),
                           overwrite = FALSE) {
  destFile <- basename(url)
  if(overwrite == TRUE) {
    download.file(url = url, destfile = destFile)
  } else if (overwrite == FALSE) {
    if(!file.exists(destFile)) {
      download.file(url = url, destfile = destFile)
    }
  }
  return(destFile)
}

#' getUniprotData
#'
#' Batch download different types of proteome-wide annotation data from Uniprot
#' using REST services
#'
#' @param outDir Target download folder path
#' @param format File format. Accepted formats are gff, fasta, txt
#' @param organism Organism code e.g. 9606 for human
#' @param reviewed Only download reviewed uniprot entries (default: TRUE)
#' @param update Boolean (default: FALSE), whether to download updated
#' annotations to the target folder
#' @return Path to the downloaded file
#' @examples
#' #download Tobacco Rattle Virus proteome (it has only 6 genes)
#' #https://www.uniprot.org/proteomes/UP000001669
#' up <- getUniprotData(outDir = getwd(), format = 'fasta',
#'                     reviewed = TRUE, organism = 652939)
#' up_fasta <- Biostrings::readAAStringSet(up)
#' \dontrun{
#' #download gff annotation for human reviewed proteome
#' gffFile <- getUniprotData(outDir = getwd(), format = 'fasta',
#'                     reviewed = TRUE, organism = 9606)
#' gff <- rtracklayer::import.gff(gffFile)
#' }
#' @export
getUniprotData <- function(outDir, format, reviewed = TRUE, organism = 9606, update = FALSE) {
  if(!format %in% c('gff', 'fasta', 'txt')) {
    stop("Error: can only download gff, fasta, or txt format files")
  }

  checkExists <- grepl(pattern = paste0('uniprot.',organism,'.*',format,'$'),
                       x = dir(outDir))

  if(sum(checkExists) > 0 & update == FALSE) {
    warning("Destination folder already contains previously downloaded annotation data: ",
            "\n See: ",dir(outDir)[checkExists],
            "\n Set 'update' to TRUE to download an updated annotation\n")
    return(dir(outDir, full.names = TRUE)[checkExists])
  }

  url <- 'https://www.uniprot.org/uniprot/?query='

  if(reviewed == TRUE) {
    url <- paste0(url, "reviewed:yes+AND+")
  }

  url <- paste0(url, "organism:",organism, "&format=",format)

  downloadDate <- paste(unlist(strsplit(date(), ' '))[c(2,4,6)], collapse = '_')
  outFile <- file.path(outDir,
                       paste(c("uniprot", organism, downloadDate, format), collapse = '.'))


  download.file(url = url, destfile = outFile)
  return(outFile)
}

#' getElmClasses
#'
#' Scrape ELM classes table from elm.eu.org
#'
#' @return A list of regular expressions from the ELM database.
#' @examples
#' elm <- getElmClasses()
#' @export
getElmClasses <- function () {
  urlData <- RCurl::getURL(url = "http://elm.eu.org/elms")
  htmlData <- XML::htmlParse(urlData)
  tables <- XML::readHTMLTable(doc = htmlData)
  elmClasses <- tables[[1]]
  colnames(elmClasses) <- gsub(pattern = ' ',
                               replacement = '_',
                               x = colnames(elmClasses))

  motifRegex <- as.list(as.vector(elmClasses$RegEx))
  names(motifRegex) <- as.vector(elmClasses$ELM_Identifier)

  return(motifRegex)
}


#' getElmInstances
#'
#' Using REST services, download ELM instances table from elm.eu.org
#' @param query (default: * -all instances-) any query term to search the
#'   database for instances
#' @param instanceLogic (default: true+positive) other options: true+negative,
#'   false+positive, false+negative
#' @param taxon (default: homo+sapiens) other option examples: mus+musculus,
#'   drosophila+melanogaster etc.
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom utils download.file
#' @importFrom utils read.table
#' @return A GRanges object where seqnames are Uniprot accessions and
#' ranges are the short linear motif locations in the proteins.
#' @examples
#' elmInst <- getElmInstances('*', 'true+positive', 'homo+sapiens')
#' @export
getElmInstances <- function (query = '*',
                             instanceLogic = 'true+positive',
                             taxon = 'homo+sapiens') {
  outFile <- file.path(getwd(), 'elmInstances.tsv')
  if(!file.exists(outFile)) {
    myUrl <- paste0('http://elm.eu.org/instances.tsv?q=',query,
                    '&instance_logic=',instanceLogic,
                    '&taxon=',taxon
                    )
    download.file(url = myUrl, destfile = outFile)
  } else {
    warning("elmInstances.tsv exists at current folder",getwd(),
            ", a new one won't be downloaded. Remove the existing
            file and re-run the function to update the file")
  }
  df <- read.table(outFile, sep = '\t', skip = 5, header = TRUE)
  elmInstances <- GenomicRanges::makeGRangesFromDataFrame(
    df = df,
    keep.extra.columns = TRUE,
    ignore.strand = TRUE,
    seqnames.field = 'Primary_Acc',
    start.field = 'Start',
    end.field = 'End')
  return(elmInstances)
}

#' getELMdomainInteractions
#'
#' Scrape known motif class - domain pairs from ELM
#' http://elm.eu.org/interactiondomains
#'
#' @importFrom RCurl getURL
#' @importFrom XML htmlParse
#' @importFrom XML readHTMLTable
#' @return A data.frame object with ELM classes versus cognate PFAM domains
#' @examples
#'
#' elm2pfam <- getELMdomainInteractions()
#'
#' @export
getELMdomainInteractions <- function () {
  urlData <- RCurl::getURL(url = "http://elm.eu.org/interactiondomains")
  htmlData <- XML::htmlParse(urlData)
  tables <- XML::readHTMLTable(doc = htmlData)
  elm2pfam <- tables[[1]]
  colnames(elm2pfam) <- gsub(pattern = ' ',
                             replacement = '_',
                             x = colnames(elm2pfam))
  return(elm2pfam[,1:3])
}

#' runIUPred
#'
#' run IUPred tool to get disorder propensity scores of a given protein sequence
#' in fasta format.
#'
#' Given a fasta file containing one or more protein sequences; calculate
#' per-base iupred disorder scores. IUPred source code can be dowloaded from
#' here: http://iupred.enzim.hu/Downloads.php After unpacking the source code,
#' cd to the src directory. Compile the code with "cc iupred.c -o iupred"
#'
#' @param iupredPath The path to the folder containing the source code of the
#'   IUPred tool
#' @param fastaFile A fasta file containing the amino acid sequences
#' @param nodeN default: 1. Positive integer for the number of parallel cores to
#'   use for downloading and processing the files
#' @importFrom parallel clusterExport
#' @importFrom Biostrings readAAStringSet
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom pbapply pblapply
#' @return A list of numerical arrays, where each value in the array correspond
#'   to the disorder score propensity of the corresponding residue in the
#'   protein sequence.
#' @examples
#' \dontrun{
#' glutFastaFile <- system.file("extdata", "glut.fasta", package = 'slimR')
#' #get disorder scores
#' iupred <- runIUPred('/home/buyar/tools/iupred/', glutFastaFile, 1)
#' #plot disorder score profile
#' plot(iupred[[1]], col = ifelse(iupred[[1]] > 0.4, 'red', 'black'), main = names(iupred)[1])
#' }
#' @export
runIUPred <- function (iupredPath,
                       fastaFile,
                       nodeN = 1) {

  fasta <- Biostrings::readAAStringSet(filepath = fastaFile, format = 'fasta')
  cl <- parallel::makeCluster(nodeN)
  parallel::clusterExport(cl = cl, varlist = c('fasta', 'iupredPath'), envir = environment())
  results <- pbapply::pblapply(cl = cl, X = 1:length(fasta), FUN = function(i) {
    requireNamespace('Biostrings')
    sequence <- fasta[i]
    #write sequence to a temp file
    id <- paste0(c('tmp_', sample(1:9, 10, replace = T)),  collapse = '')
    tmpFasta <- file.path(getwd(), paste0(id, '.fasta'))
    Biostrings::writeXStringSet(x = sequence, filepath = tmpFasta, format = 'fasta')

    #run iupred
    tmpOut <- file.path(getwd(), paste0(id, '.iupred'))
    myCommand <- paste(paste0("export IUPred_PATH=",iupredPath, ";"),
                       file.path(iupredPath, 'iupred'),
                       tmpFasta, "short >", tmpOut)

    system(command = myCommand)

    if(file.exists(tmpOut)) {
      df <- read.table(file = tmpOut)
      result <- as.numeric(df$V3)
    } else {
      result <- NA
    }
    #clean up tmp files
    file.remove(tmpFasta)
    file.remove(tmpOut)
    return(result)
  })
  parallel::stopCluster(cl)
  names(results) <- names(fasta)
  return(results)
}


#' getPFAM
#'
#' Download pfam domain annotations of complete proteomes
#'
#' @param organism Organism taxonomy id (e.g. human: 9606)
#' @param pfam_version Version of PFAM release to use: default Pfam30.0. Check
#'   ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/ for available versions.
#' @importFrom data.table fread
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom utils download.file
#' @return A Granges object containing the coordinates PFAM domains in protein
#'   sequences
#' @examples
#' \dontrun{
#' pfam <- getPFAM(organism = 9606, pfam_version = 'Pfam30.0')
#' }
#' @export
getPFAM <- function(organism = 9606, pfam_version = "Pfam30.0") {

  url <- paste0('ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/',
                pfam_version,
                '/proteomes/',
                organism,
                ".tsv.gz")

  outFile <- file.path(getwd(), paste0(pfam_version, '.', organism, '.tsv.gz'))

  if(file.exists(outFile)) {
    warning("PFAM annotation data already exists at: \n \t",outFile,"\n")
  } else {
    download.file(url = url, destfile = outFile)
  }

  pfam <- data.table::fread(outFile, sep = '\t', header = F)
  colnames(pfam) <- c('seqname', 'start', 'end', 'env_start', 'env_end',
                      'pfam_acc', 'pfam_name', 'pfam_type',
                      'hmm_start', 'hmm_end', 'hmm_length', 'bit_score',
                      'e_value', 'clan'
                      )
  pfam <- GenomicRanges::makeGRangesFromDataFrame(pfam,
                                                  keep.extra.columns = TRUE)

  return(pfam)
}

#' getPFAMClans
#'
#' Download PFAM domain clans
#'
#' @importFrom RCurl getURL
#' @importFrom XML htmlParse
#' @importFrom XML readHTMLTable
#' @return A data.frame object
#' @examples
#' \dontrun{
#' clans <- getPFAMClans()
#' }
#' @export
getPFAMClans <- function() {
  urlData <- RCurl::getURL(url = "http://pfam.xfam.org/clan/browse")
  htmlData <- XML::htmlParse(urlData)
  tables <- XML::readHTMLTable(doc = htmlData)
  clans <- tables$clanBrowse
  return(clans)
}


#' validateVariants
#'
#' This function is used to validate whether the missense variants are consistent
#' with the fasta sequences of the proteins
#'
#' @param df A data.frame or data.table containing minimally the following
#'   columns: 1. uniprotAccession 2. wtAA (wild-type amino acid) 3. pos (the position
#'   of the wild-type amino acid in the protein sequence)
#' @param fasta AAStringSet object of protein sequences, in which names of the list items
#'   should correspond to the uniprotAccession column of the 'df' input
#' @param nodeN Number of cores to use to parallelise the run
#' @return A data.frame consisting of the input data.frame and two additional columns
#' for
#' 1. mappedResidue (the actual residue at the variant position)
#' 2. validity (boolean TRUE/FALSE) showing if the mapped residue matches the
#' wild-type residue in the corresponding variant.
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel stopCluster
#' @importFrom pbapply pbapply
#' @examples
#' glutFastaFile <- system.file("extdata", "glut.fasta", package = 'slimR')
#' glutFasta <- Biostrings::readAAStringSet(glutFastaFile)
#' data("glutMutations")
#' val <- validateVariants(df = glutMutations[1:10,], fasta = glutFasta, nodeN = 1)
#' @export
validateVariants <- function(df, fasta, nodeN = 1) {
  cl <- parallel::makeCluster(nodeN)
  parallel::clusterExport(cl, varlist = c('df', 'fasta'),  envir = environment())
  mappedResidues <- pbapply::pbapply(cl = cl, X = df, MARGIN = 1, FUN = function(x) {
    requireNamespace('Biostrings')
    uni <- x[['uniprotAccession']]
    pos <- as.numeric(x[['pos']])
    ifelse(pos > length(fasta[[uni]]), NA, as.character(fasta[[uni]][pos]))
  })
  parallel::stopCluster(cl)

  df$mappedResidue <- mappedResidues
  df$validity <- df$wtAA == df$mappedResidue
  return(df)
}

