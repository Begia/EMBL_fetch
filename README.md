# EMBL_fetch
Fetch protein and nucleotide records from EMBL databases (work in progress)

Since getting the actual DNA sequence from NCBI and ENA archives can sometimes be difficult, I made a tool
that extracts the required region of the nucleotide record, protein record and some useful metadata.

The program can be run as a standalone R script. User provides an input file with protein and/or nucleotide 
accession numbers (single column list). The program then goes through the list and queries various EMBL 
databases for records of the required accession numbers, finds the original nucleotide record and extracts 
the region required. Also upstream and downstream options are available. The program has only been tested 
with bacterial records. Introns in eukaryotic DNA might cause issues.
