# Short Linear Motif (SLiM) Analysis in the context of human diseases

slimR is a protein sequence analysis package centered around
  short linear motifs (SLiMs) and their connections to human diseases.
  The package contains functions to retrieve data from public databases
  such as ELM (elm.eu.org), UniProt (uniprot.org). Annotated SLiMs and
  SLiM patterns  (regular expressions) are retrieved from the ELM database.
  Protein sequence  features such as disease-causing mutations, polymorphisms,
  protein domains,  and many other annotated features are retrieved from the
  UniProt database.  The main motivation of the package is to facilitate
  analysis of the impact of protein sequence variations due to mutations
  and polymorphism on the content of SLiMs. Besides functions to retrieve
  data from public resources, there are functions to calculate disorder
  scores, search for SLiM instances in given sequences, look for gain/loss
  of SLiMs due to sequence variations. The package should be useful not
  only for the analysis of human diseases but also studies of the evolution
  of species via gain/loss of SLiMs in any comparative context.

## External dependency 
Currently, the package depends on IUPred tool for sequence disorder score prediction.
IUPred source code can be dowloaded from here: http://iupred.enzim.hu/Downloads.php 
After unpacking the source code, change to the src directory. 
Compile the code with "cc iupred.c -o iupred". 
Some of the functions will require the path to the folder that contains the iupred executable

## Usage

`createDB(uniprotAccessions = c('P04637', 'P11166', 'P06400'))`

Given a vector of uniProt accessions, extract and process all available
sequence features and save the resulting objects in folders and RDS objects
for: 

1. Fasta sequences 
2. Disorder predictions using IUPred 
3. Uniprot feature files downloaded in gff format 
4. Uniprot variants including disease-causing and polymorphic substitutions 
5. Short linear motifs - all substrings matching all available patterns 
6. SLiM Changes due to disease causing mutations - The collection of changes in the proteins' SLiM content when diesase-causing variants at 4) are applied to the sequence 
7. SLiM Changes due to polymorphisms - The collection of changes in the proteins' SLiM content when polymorphisms at 4) are applied to the sequence

The function can be run using multiple cores using the `nodeN` parameter. This would be especially useful if a database is created for the whole proteome (human proteome contains >20000 sequences)
for the first time. Once the database is created, it will be very quick to do integrated queries. 




