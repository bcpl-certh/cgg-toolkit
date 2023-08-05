# MagicMatch--cross-referencing sequence identifiers across databases

original authors: Smith M, Kunin V, Goldovsky L, Enright AJ, Ouzounis CA •
PMID: 15961438 DOI: 10.1093/bioinformatics/bti548

updates by: Dimitris Vasiliou, Christos A Ouzounis, [BCPLab](http://genome.academy/) 2023


## Motivation: 
At present, mapping of sequence identifiers across databases is a daunting, time-consuming and computationally expensive process, usually achieved by sequence similarity searches with strict threshold values.

## Summary:
We present a rapid and efficient method to map sequence identifiers across databases. The method uses the MD5 checksum algorithm for message integrity to generate sequence fingerprints and uses these fingerprints as hash strings to map sequences across databases. The program, called MagicMatch, is able to cross-link any of the major sequence databases within a few seconds on a modest desktop computer.

### Build (Linux)

Make sure that gcc and g++ are installed.

```
make
./magicmatch 
```

### Usage: 

  ```
  ./magicmatch [-help show usage] [-f file1 file2 ...] [-fe file1 file2 ...] [-t translate to ID-KEY || -c find identicals between files]
  
  where:
  	-f	fasta file names (at least one file is required) 
  	-fe	previously encoded fasta files 
  	-t	translate the sequences of the files 
  	-c	map between ids of all the files
  	-h	show usage
  ```
