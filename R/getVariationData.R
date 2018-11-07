

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
#' @return A subset of the initial data.frame or data.table object such that the
#'   given positions of amino acids match the residues in the fasta sequences
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel stopCluster
#' @export
validateVariants <- function(df, fasta, nodeN = 8) {
  cl <- parallel::makeCluster(10)
  parallel::clusterExport(cl, varlist = c('fasta'))
  mappedResidues <- pbapply::pbapply(cl = cl, X = df, MARGIN = 1, FUN = function(x) {
    require(Biostrings)
    uni <- x[['uniprotAccession']]
    pos <- as.numeric(x[['pos']])
    ifelse(pos > length(fasta[[uni]]), NA, as.character(fasta[[uni]][pos]))
  })
  parallel::stopCluster(cl)

  df$mappedResidue <- mappedResidues
  df$validity <- df$wtAA == df$mappedResidue
  return(df)
}


