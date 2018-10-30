

#' runVEP
#'
#' A wrapper function to run Ensembl's variant_effect_predictor script
#'
#' @param perlPATH path to PERL installation to use
#' @param vepPATH path to the variant_effect_predictor.pl script
#' @param inputFileName path to input file that contains variation data in a
#'   format acceptable to variant_effect_predictor software (see:
#'   http://www.ensembl.org/info/docs/tools/vep/vep_formats.html#default)
#' @param outputFileName file name to write the results
#' @param overwrite (default: FALSE) set it to TRUE to overwrite the existing
#'   VEP output file.
#' @return a data.table data.frame containing variation data read from VEP
#'   output
#'
#' @importFrom data.table fread
#' @export
runVEP <- function(vepPATH = '/home/buyar/.local/bin/variant_effect_predictor.pl',
                   inputFileName,
                   outputFileName = 'VEPresults.tsv',
                   overwrite = FALSE,
                   nodeN = 4) {
  if(!file.exists(outputFileName) | (file.exists(outputFileName) & overwrite == TRUE)) {
    command <- paste(vepPATH,'-i',inputFileName,' -o',outputFileName,
                     ' --cache --uniprot --force_overwrite --fork',nodeN)
    if(!is.null(perlLibPath)) {
      command <- paste0("export PERL5LIB=",perlLibPath,";",command)
    }
    message("Running VEP with the command:\n",command)
    system(command)
  }
  return(outputFileName)
}

#' processVEP
#'
#' This function processes the output of Variant Effect Predictor to select
#' missense variants and create some columns that are useful to assess the
#' pathogenicity of variants
#'
#' @param vepFilePath path to the VEP results obtained from running
#'   variant_effect_predictor on the given vcfFilePath
#' @param nodeN (default: 4) Number of cores to use for parallel processing
#' @importFrom data.table fread
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel stopCluster
#' @return A data.table object
#' @export
processVEP <- function(vepFilePath, overwriteRDS = FALSE, nodeN = 4) {
  if(!file.exists(vepFilePath)) {
    stop("Couldn't find the path to the VEP results file ",vepFilePath)
  }

  #read VEP results
  vepRDS <- paste0(vepFilePath, '.RDS')
  if(file.exists(vepRDS) & overwriteRDS == FALSE) {
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
    colnames(vep)[1] <- 'Identifier'

    #remove rows where multiple amino acids are reported for missense variants
    vep <- vep[grep('^.\\/.$', vep$Amino_acids),]
    cat("Removed variants where multiple amino acids are reported for missense variants...\n",
        "\tkeeping",nrow(vep),"mappings of",nrow(unique(vep[,1])),
        "unique variants from",vepFilePath,"\n")

    #add extra columns about the mutation positions and amino acids
    vep$pos <- as.numeric(vep$Protein_position)
    vep$wtAA <- gsub(pattern = '\\/.$', '', vep$Amino_acids)
    vep$mutAA <- gsub(pattern = '^.\\/', '', vep$Amino_acids)

    saveRDS(object = vep, file = vepRDS)
  }
  return(vep)
}

#' getClinVarVEPmissenseVariants
#'
#' This function downloads ClinVar variants and runs variant_effect_predictor
#' (VEP) tool and merges the clinvar variants with missense variants predicted
#' by VEP.
#'
#' @return A data.table object
#' @export
getClinVarVEPmissenseVariants <- function(vepPath = '/home/buyar/.local/bin/variant_effect_predictor.pl',
                                          clinvarDataURL, overwrite = FALSE, nodeN = 4) {

  rdsFile <- file.path(getwd(), 'clinvar_VEP_missenseVariants.RDS')
  if(file.exists(rdsFile) & overwrite == FALSE) {
    return(readRDS(rdsFile))
  } else {
    clinvarDataFile <- getClinVarData(url = clinvarDataURL,
                                      overwrite = overwrite)
    clinvarData <- data.table::fread(clinvarDataFile)

    clinvarData <- clinvarData[Type == 'single nucleotide variant' &
                                 Assembly == 'GRCh38']

    #subset the clinvar data by columns, write to a file and pass to runVEP
    #function
    clinvarData$Allele <- paste(clinvarData$ReferenceAllele,
                                clinvarData$AlternateAllele, sep = '/')
    #clinvar uses the positive strand as reference nucleotide
    clinvarData$Strand <- '+'
    clinvarData$Name <- gsub(' ', '_', clinvarData$Name)
    # identifier created to uniquely refer to every row in the data table. this
    # is used as input to VEP and later used to merge VEP results with some of
    # the clinvarData fields (e.g. dbSNP id etc)
    clinvarData$Identifier <- paste(clinvarData$Chromosome, clinvarData$Start,
                                    clinvarData$Stop, clinvarData$Name,
                                    sep = ':')
    write.table(x = clinvarData[,c('Chromosome', 'Start', 'Stop',
                                   'Allele', 'Strand', 'Identifier')],
                file = 'clinvarData.processed.tsv', quote = F,
                sep = '\t', row.names = F, col.names = F)

    runVEP(vepPATH = vepPath,
           inputFileName = 'clinvarData.processed.tsv',
           outputFileName = 'clinvarData.processed.VEPoutput.tsv',
           overwrite = TRUE,
           nodeN = nodeN)

    vep <- processVEP(vepFilePath = 'clinvarData.processed.VEPoutput.tsv',
                      overwriteRDS = TRUE,
                      nodeN = nodeN)

    #merge vep data with clinvar data (notice that vep results will only contain
    #missense-variants, so the clinvar data will be down-sized to only those
    #that have a calculated consequence of 'missense variant' according to vep)
    clinvarVEPdata <- merge(clinvarData, vep, by = 'Identifier')
    saveRDS(object = clinvarVEPdata,
            file = rdsFile)
    return(clinvarVEPdata)
  }
}


#' combineClinVarWithHumsavar
#'
#' This function processes humsavar variants (output of getHumSavar()) and
#' clinvar variants (output of runVEP) and merges into a simplified data.table
#' object
#'
#' @param clinvarVEPdata path to VCF file containing variation data from ClinVar
#'   database
#' @return A data.table object
#'
#' @importFrom data.table data.table
#' @export
combineClinVarWithHumsavar <- function(clinvarVEPdata) {
  #1. reformat clinvar data
  cv <- unique(subset(clinvarVEPdata,
                      select = c('uniprotAccession', 'pos', 'ClinicalSignificance',
                                 'RS# (dbSNP)', 'wtAA', 'mutAA', 'PhenotypeList')))
  colnames(cv) <- c('uniprotAccession', 'pos', 'variant', 'dbSNP',
                    'wtAA', 'mutAA', 'clinvarDisease')

  #further process and filter clinvar data - only keep entries with good
  #confidence of being either disease-related or confidently benign
  cv <- cv[unlist(lapply(cv$variant, function(x) sum(c('Pathogenic', 'Pathogenic/Likely pathogenic', 'Benign') %in% unlist(strsplit(x, ','))) > 0))]
  #remove pathogenic variants that don't have a specified disease name
  cv$clinvarDisease <- unlist(lapply(cv$clinvarDisease,
                                     function(x) {
                                       update <- setdiff(unlist(strsplit(x, ';')),
                                                         c('not specified', 'not provided'))
                                       if(length(update) > 1){
                                         return(paste0(update, collapse = '; '))
                                       } else if (length(update) == 1){
                                         return(update)
                                       } else {
                                         return(NA)
                                       }}))
  cv <- cv[!(grepl('Pathogenic', cv$variant, ignore.case = T) & is.na(clinvarDisease))]

  cv$dbSNP <- ifelse(cv$dbSNP == -1, '-', paste0('rs', cv$dbSNP))

  #map terms to 'Disease' or 'Polymorphism'
  cv[grepl('Pathogenic', cv$variant, ignore.case = T)]$variant <- 'Disease'
  cv[grepl('Benign', cv$variant, ignore.case = T)]$variant <- 'Polymorphism'

  #2. get humsavar data
  hs <- data.table::data.table(as.data.frame(getHumSavar()))
  hs$variant <- as.character(hs$variant)
  #humsavar data simplified
  hs <- unique(subset(hs, select = c('seqnames', 'start', 'variant',
                                     'dbSNP', 'wtAA', 'mutAA', 'diseaseName')))
  colnames(hs) <- c('uniprotAccession', 'pos', 'variant', 'dbSNP',
                    'wtAA', 'mutAA', 'humsavarDisease')

  #subset to get only disease or polymorphism variants
  hs <- hs[variant %in% c('Disease', 'Polymorphism')]

  combined <- merge(hs, cv, by = c('dbSNP', 'uniprotAccession', 'pos', 'wtAA', 'mutAA'), all = T)
  colnames(combined) <- c('dbSNP', 'uniprotAccession', 'pos', 'wtAA', 'mutAA',
                          'humsavarVariant', 'humsavarDisease',
                          'clinvarVariant', 'clinvarDisease')

  combined <- unique(combined)

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
          flankN <- residues[(pos-flankLength):(pos-1)]
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


