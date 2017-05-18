#' parseMutation
#'
#' Given a vector of mutation substitutions (e.g. "p.His160Arg")
#' -> split "p.His160Arg" into "H 160 R"
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
#' @return A Granges object containing the coordinates of mutated sites in proteins
#' @export
getHumSavar <- function () {
  variantFile <- file.path(getwd(), 'humsavar.txt')
  if (!file.exists(variantFile)) {
    download.file(url = 'www.uniprot.org/docs/humsavar.txt',
                  destfile = variantFile, method = "curl")
  } else {
    warning("humsavar.txt exists at current folder",getwd(),
            ", a new one won't be downloaded. Remove the existing
            file and re-run the function to update the file")
  }

  if(file.exists(paste0(variantFile, '.RDS'))){
    return(readRDS(paste0(variantFile, '.RDS')))
  } else {
    #skip first 30 lines which don't contain mutation data
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
    saveRDS(object = mut, file = paste0(variantFile, ".RDS"))
    return(mut)
  }
  }

#' parseUniprotHumanVariationData
#'
#' Parse the human variation data (homo_sapiens_variation.txt.gz) from Uniprot.
#' The variation annotation contains polymorphism data from COSMIC, 1000GP etc.
#' and mapped to Uniprot proteins.
#' See ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/variants/homo_sapiens_variation.txt.gz
#' TODO: add a parameter to overwrite the rds file generated
#' @param filePath The path to the human variation data file (e.g homo_sapiens_variation.txt.gz) downloaded from Uniprot
#' @return A data.table format table of mutation data
#' @export
parseUniprotHumanVariationData <- function (filePath, outFile = 'parseUniprotHumanVariationData.tsv',
                                            keepColumns = c(2,3,5,6,7,8,14),
                                            sourceDB = '1000Genomes',
                                            consequenceType = 'missense variant') {
  if (!file.exists(paste0(filePath, '.rds'))) {
    con <- file(filePath, 'r')
    out <- file(outFile, 'w')
    while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
      fields <- unlist(strsplit(oneLine, '\t'))
      if(length(fields) == 14) {
        if (fields[5] == consequenceType & grepl(sourceDB, fields[14]) == TRUE) {
          writeLines(text = paste(fields[keepColumns], collapse = '\t'), con = out)
        }
      }
    }
    close(con = con)
    close(con = out)
    dt <- data.table::fread(outFile, sep = '\t', header = FALSE)
    saveRDS(object = dt, file = paste0(filePath, '.rds'))
  } else {
    dt <- readRDS(paste0(filePath, '.rds'))
  }
  return(dt)
}

### Functions to download and map ClinVar data to Uniprot sequences
#' getClinVarData
#'
#' This function will fetch the clinvar data from the given url and parse the contents
#' of the downloaded file.
#'
#' @param url The url to the ftp location of the clinvar dataset
#'  (e.g. ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20170404.vcf.gz)
#' @return A tibble object extracted from the downloaded vcf file
#' @importFrom vcfR read.vcfR
#' @export
getClinVarData <- function(url, overwrite = FALSE, parseDownloadedFile = TRUE) {
  destFile <- basename(url)
  if(overwrite == TRUE) {
    download.file(url = url, destfile = destFile)
  } else if (overwrite == FALSE) {
    if(!file.exists(destFile)) {
      download.file(url = url, destfile = destFile)
    }
  }

  if(file.exists(destFile)) {
    if(overwrite == TRUE) {
      gunzipCommand <- paste('gunzip -f',destFile)
    } else {
      gunzipCommand <- paste('gunzip',destFile)
    }
    system(gunzipCommand)
    destFile <- gsub('.gz$', '', destFile)
    if(parseDownloadedFile == TRUE) {
      clinvarData <- vcfR::extract_info_tidy(vcfR::read.vcfR(destFile))
      return(clinvarData)
    }
  } else {
    stop("Couldn't find",destFile,"to parse the results.
         Probably the download didn't work\n")
  }
}

#' runVEP
#'
#' A wrapper function to run Ensembl's variant_effect_predictor script for a
#' given vcf file.
#'
#' @param vepPATH path to the variant_effect_predictor.pl script
#' @param vcfFilePath path to VCF file containing variation data
#'
#' @return a data.table data.frame containing variation data read from VEP
#'   output
#'
#' @importFrom data.table fread
#' @export
runVEP <- function(vepPATH = '/home/buyar/.local/bin/variant_effect_predictor.pl',
                   vcfFilePath,
                   overwrite = FALSE,
                   parseResultFile = FALSE) {
  vepOutFile <- gsub(pattern = '.vcf$', replacement = '.VEP.tsv', x = vcfFilePath)
  if(!file.exists(vepOutFile) | (file.exists(vepOutFile) & overwrite == TRUE)) {
    system(paste(vepPATH,'-i',vcfFilePath,' -o',vepOutFile,' --cache --uniprot --force_overwrite'))
  }
  if(parseResultFile == TRUE) {
    vepData <- data.table::fread(vepOutFile)
    return(vepData)
  }
}

#' processVEP
#'
#' This function processes the output of Variant Effect Predictor to select
#' missense variants and create some columns that are useful to assess the
#' pathogenicity of variants
#'
#' @param vcfFilePath path to VCF file containing variation data
#' @param vepFilePath path to the VEP results obtained from running
#'   variant_effect_predictor on the given vcfFilePath
#' @param nodeN (default: 8) Number of cores to use for parallel processing
#' @importFrom data.table fread
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel stopCluster
#' @importFrom vcfR read.vcfR
#' @importFrom vcfR extract_info_tidy
#' @return A data.table object
#' @export
processVEP <- function(vcfFilePath, vepFilePath, nodeN = 4) {
  if(!file.exists(vcfFilePath)) {
    stop("Couldn't find the path to the vcf file ",vcfFilePath)
  }

  if(!file.exists(vepFilePath)) {
    stop("Couldn't find the path to the VEP results file ",vepFilePath)
  }

  #read clinvar data
  clinvarRDS <- gsub(pattern = '.vcf', replacement = '.RDS', x = vcfFilePath)
  if(file.exists(clinvarRDS)) {
    clinvarData <- readRDS(clinvarRDS)
  } else {
    clinvarData <- vcfR::extract_info_tidy(vcfR::read.vcfR(vcfFilePath))
    saveRDS(object = clinvarData, file = clinvarRDS)
  }

  #read VEP results
  vepRDS <- gsub(pattern = '.tsv', replacement = '.RDS', x = vepFilePath)
  if(file.exists(vepRDS)) {
    vep <- readRDS(vepRDS)
  } else {
    vepRaw <- data.table::fread(vepFilePath)
    cat("Read",nrow(vepRaw),"mappings of",nrow(unique(vepRaw[,1])),"unique variants from",vepFilePath,"\n")

    #filter for missense variants with swissprot ids
    vep <- vepRaw[grepl(pattern = 'SWISSPROT', x = vepRaw$Extra)]

    cat("Removing mappings that don't have a SWISSPROT ID...\n",
        "\tkeeping",nrow(vep),"mappings of",nrow(unique(vep[,1])),
        "unique variants from",vepFilePath,"\n")

    vep <- vep[Consequence == 'missense_variant']

    cat("Removing variants that don't lead to a single amino-acid substitution...\n",
        "\tkeeping",nrow(vep),"mappings of",nrow(unique(vep[,1])),
        "unique variants from",vepFilePath,"\n")

    cl <- parallel::makeCluster(nodeN)
    parallel::clusterExport(cl, varlist = c('vep'), envir = environment())
    vep$uniprotAccession <- gsub(pattern = '(SWISSPROT=|;$)', replacement = '', do.call(c, parLapply(cl, vep$Extra, function(x) {
      unlist(stringi::stri_extract_all(str = x, regex = 'SWISSPROT=.*?;'))
    })))
    parallel::stopCluster(cl)
    colnames(vep)[1] <- 'dbSNP'

    #remove rows where multiple amino acids are reported for missense variants
    vep <- vep[grep('^.\\/.$', vep$Amino_acids),]
    cat("Removed variants where multiple amino acids are reported for missense variants...\n",
        "\tkeeping",nrow(vep),"mappings of",nrow(unique(vep[,1])),
        "unique variants from",vepFilePath,"\n")

    #add extra columns about the mutation positions and amino acids
    vep$pos <- as.numeric(vep$Protein_position)
    vep$wtAA <- gsub(pattern = '\\/.$', '', vep$Amino_acids)
    vep$mutAA <- gsub(pattern = '^.\\/', '', vep$Amino_acids)
    vep$RS <- as.numeric(gsub('rs', '', vep$dbSNP))

    #add extra info from clinvar data
    vep$CLNSIG <- clinvarData[match(vep$RS, clinvarData$RS),]$CLNSIG
    vep$CLNDBN <- clinvarData[match(vep$RS, clinvarData$RS),]$CLNDBN
    vep$isCommonVariant <- clinvarData[match(vep$RS, clinvarData$RS),]$COMMON
    vep$pathogenic <- unlist(lapply(vep$CLNSIG, function(x) sum(c('4', '5') %in% unlist(strsplit(x,  split = '\\|'))) > 0))
    vep$implicatedInAnyDisease <- gsub(pattern = "not_provided|not_specified|\\|", replacement = '', vep$CLNDBN) != ''

    saveRDS(object = vep, file = vepRDS)
  }
  return(vep)
}

#' combineClinVarWithHumsavar
#'
#' This function processes humsavar variants (output of getHumSavar()) and
#' clinvar variants (output of runVEP) and merges into a simplified data.table
#' object
#'
#' @param vcfFilePath path to VCF file containing variation data from ClinVar
#'   database
#' @param vepFilePath path to the VEP results obtained from running
#'   variant_effect_predictor on the given vcfFilePath
#' @param nodeN (default: 8) number of cpu nodes to use when running in parallel
#' @return A data.table object
#'
#' @importFrom data.table data.table
#' @export
combineClinVarWithHumsavar <- function(vcfFilePath, vepFilePath, nodeN = 8) {
  #clinvar VEP results simplified
  cv <- processVEP(vcfFilePath, vepFilePath, nodeN = nodeN)
  cv$variant <- 'Unclassified'
  cv[cv$pathogenic == TRUE]$variant <- 'Disease'
  cv[cv$pathogenic == FALSE & cv$isCommonVariant == 1,]$variant <- 'Polymorphism'

  cv <- unique(subset(cv, select = c('uniprotAccession', 'pos', 'variant', 'dbSNP', 'wtAA', 'mutAA', 'CLNDBN')))
  colnames(cv) <- c('uniprotAccession', 'pos', 'variant', 'dbSNP', 'wtAA', 'mutAA', 'clinvarDisease')

  hs <- data.table::data.table(as.data.frame(getHumSavar()))

  #humsavar data simplified
  hs <- unique(subset(hs, select = c('seqnames', 'start', 'variant', 'dbSNP', 'wtAA', 'mutAA', 'diseaseName')))
  colnames(hs) <- c('uniprotAccession', 'pos', 'variant', 'dbSNP', 'wtAA', 'mutAA', 'humsavarDisease')

  combined <- merge(hs, cv, by = c('dbSNP', 'uniprotAccession', 'pos', 'wtAA', 'mutAA'), all = T)
  colnames(combined) <- c('dbSNP', 'uniprotAccession', 'pos', 'wtAA', 'mutAA',
                          'humsavarVariant', 'humsavarDisease', 'clinvarVariant', 'clinvarDisease')

  return(combined)
}

#' validateVariants
#'
#' This function is used to validate whether the missense variants are consistent
#' with the fasta sequences of the proteins
#'
#' @param df A data.frame or data.table containing minimally the following
#'   columns: 1. uniprotAccession 2. wtAA (wild-type amino acid) 3. pos (the position
#'   of the wild-type amino acid in the protein sequence)
#' @param fasta A list of fasta sequences, in which names of the list items
#'   should correspond to the uniprotAccession column of the 'df' input
#' @param nodeN Number of cores to use to parallelise the run
#' @return A subset of the initial data.frame or data.table object such that the
#'   given positions of amino acids match the residues in the fasta sequences
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel stopCluster
#' @export
validateVariants <- function(df, fasta, nodeN = 8) {
  flankLength <- 7
  cl <- parallel::makeCluster(nodeN)
  parallel::clusterExport(cl, varlist = c('df', 'fasta', 'flankLength'),
                          envir = environment())
  validity <- data.frame(do.call(rbind, parLapply(cl, 1:nrow(df), function(i) {
    uni <- df$uniprotAccession[i]
    AA <- df$wtAA[i]
    pos <- as.numeric(df$pos[i])
    status <- ''
    flankingSequence <- ''
    if(!uni %in% names(fasta)) {
      status <- 'uni_not_available_as_fasta'
    } else {
      residues <- unlist(strsplit(fasta[[uni]], ''))
      if(length(residues) < pos) {
        status <- 'invalid'
      } else if (residues[pos] == AA) {
        status <- 'valid'
        if(pos > flankLength) {
          flankN <- residues[c((pos-flankLength):pos)]
        } else {
          flankN <- residues[1:(pos-1)]
        }

        if(length(residues) - flankLength >= pos) {
          flankC <- residues[(pos+1):(pos+flankLength)]
        } else {
          flankC <- residues[(pos+1):length(residues)]
        }
        flankingSequence <- paste0(c(flankN, AA, flankC), collapse = '')
      } else {
        status <- 'invalid'
      }
    }
    return(c(status, flankingSequence))
  })), stringsAsFactors = F)
  colnames(validity) <- c('validity', 'flankingSequence')
  df <- cbind(df, validity)
  stopCluster(cl)
  return(df)
}


