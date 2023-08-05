# cogent_utils--utilities to create CoGenT-style identifier formats for computational genomics

_note_: CoGenT/CoGenT++ is a data environment for computational research in comparative and functional genomics, designed to address issues of consistency, reproducibility, scalability and accessibility; the original data collections are not provided, just some utilities that allow creation of custom databases.

original authors: Goldovsky L, Janssen P, Ahrén D, Audit A, Cases I, Darzentas N, J Enright AJ, López-Bigas N, Peregrin-Alvarez JM, Smith M, Tsoka S, Kunin V, Ouzounis CA • PMID: 16216832 DOI: 10.1093/bioinformatics/bti579

original authors: Janssen P, Enright AJ, Audit B, Cases I, Goldovsky L, Harte N, Kunin V, Ouzounis CA • PMID: 12874064 DOI: 10.1093/bioinformatics/btg161

utilities by: Dimitris Vasiliou, Christos A Ouzounis, [BCPLab](http://genome.academy/) 2023

### create_cogent.sh

create_cogent.sh is a shell script that takes as input a text file
that contains the names of protein fasta files and given titles in columns, e.g.
   
    GCF_008822105.2_bTaeGut2.pat.W.v2_protein.faa   Taeg-2p1
    GCF_015227675.2_mRatBN7.2_protein.faa   Ratn-n72
    GCF_900215245.1_IMG-taxon_2617270901_annotated_assembly_protein.faa Taxo-un3
  
and a destination folder.

When executed, it outputs the code for a new shell script.
That new script, using awk and sed magic, creates the destination folder, copies the fasta files mentioned in 
the input file to the destination folder and transforms the headers of the protein sequences in the fasta files
to a more intuitive form that includes the new title (2nd column), concatenated with an increasing counter and the previous 
protein sequence header, e.g.

Before: First two sequences headers in file GCF_008822105.2_bTaeGut2.pat.W.v2_protein.faa:
  
    >NP_001041718.1 alpha-synuclein [Taeniopygia guttata]
    >NP_001041719.1 neurocalcin-delta [Taeniopygia guttata]

After: The first two sequence headers in the generated file Taeg-2p1.faa in the destination folder:

    >Taeg-2p1-01-000000 NP_001041718.1 alpha-synuclein [Taeniopygia guttata]
    >Taeg-2p1-01-000001 NP_001041719.1 neurocalcin-delta [Taeniopygia guttata]
    
It is assumed that the current directory includes the input text file and the input fasta files mentioned in it.
    
#### Usage:
```
Usage: create_cogent.sh arg1 arg2

Creates a cogent generator script

Available options:

arg1            Input file that contains the fasta files to be transformed
arg2            Name of the generated script's output folder
-h, --help      Print this help and exit
```

```
 ./create_cogent.sh destination_folder > gen.sh
 chmod +x gen.sh
 ./gen.sh
 ```
That creates the destination folder with the new fasta files with the transformed headers.

