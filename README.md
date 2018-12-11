# Short Linear Motif (SLiM) Analysis in the context of human diseases

slimR is a protein sequence analysis package centered around short linear motifs
(SLiMs) and their connections to human diseases. The package contains functions
to retrieve data from public databases such as ELM (elm.eu.org), UniProt
(uniprot.org), PFAM, and Clinvar. Annotated SLiMs and SLiM patterns  (regular
expressions) are retrieved from the ELM database. Protein sequence  features
such as disease-causing mutations, polymorphisms, protein domains,  and many
other annotated features are retrieved from the UniProt database.  The main
motivation of the package is to facilitate analysis of the impact of protein
sequence variations due to mutations and polymorphism on the content of SLiMs.
Besides functions to retrieve data from public resources, there are functions to
calculate disorder scores, search for SLiM instances in given sequences, look
for gain/loss of SLiMs due to sequence variations. The package should be useful
not only for the analysis of human diseases but also studies of the evolution of
species via gain/loss of SLiMs in any comparative context.

# Installation

## from github
The package can be installed via github using R library `devtools`. 

> install.packages('devtools')

> devtools::install_github('BIMSBbioinfo/slimR')

## using Conda

The following command installs the latest version of the R package from Anaconda. 
> conda install -c bora.uyar r-slimr

## External dependency 
Currently, the package depends on IUPred tool for sequence disorder score prediction.
IUPred source code can be dowloaded from here: http://iupred.enzim.hu/Downloads.php 
After unpacking the source code, change to the src directory. 
Compile the code with "cc iupred.c -o iupred". 
Some of the functions will require the path to the folder that contains the iupred executable



## Acknowledgements

- [The Bioinformatics Platform](http://bioinformatics.mdc-berlin.de)

- [Berlin Institute for Medical Systems Biology](https://www.mdc-berlin.de/13800178/en/bimsb)

- [Max-Delbrueck-Center for Molecular Medicine](https://www.mdc-berlin.de)

- [The German Network for Bioinformatics Infrastructure](http://www.denbi.de/)
