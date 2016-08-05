#' Amino Acid table
#'
#' Amino acid table for converting between three letter code and one letter code.
#'
#' @return A data.frame object
#' @format data.frame object with 22 rows and two columns
"aaTable"


#' Sample sequence extracted from a fasta file
#'
#' The amino acid sequence for P11166 protein (GLUT1)
#'
#' @return A data.frame object
#' @source Downloaded from uniprot \url{www.uniprot.org/uniprot/P11166.fasta}
#' @format data.frame object with 22 rows and two columns
"sampleFasta"

#' GLUT1 disease variants
#'
#' Table of disease associated substitutions parsed from Humsavar (UniProt)
#' @return A data.frame object
#' @source GLUT1 protein \url{www.uniprot.org/uniprot/P11166}
#' and humsavar \url{www.uniprot.org/docs/humsavar.txt}
"sampleVariants"

#' SLiM regular expressions from ELM database
#'
#' List of available SLiM classes and their regular expression defined in the
#' ELM database as of 5th of August, 2016
#'
#' This list can be updated using the getElmClasses() function motifRegex <-
#' as.list(as.vector(elms$RegEx)) names(motifRegex) <-
#' as.vector(elms$ELM_Identifier)
#'
#' @return A list object
#' @source ELM classes \url{http://elm.eu.org/elms}
"motifRegex"
