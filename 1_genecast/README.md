# genecast--compositional bias filtering algorithm

original authors: Promponas V, Enright AJ, Tsoka S, Kreil DP, Leroy C, Hamodrakas S, Sander C, Ouzounis CA •
PMID: 11120681 DOI: 10.1093/bioinformatics/16.10.915

updates by: Dimitris Vasiliou, Christos A Ouzounis, [BCPLab](http://genome.academy/) 2023


 #### Usage:
 ```
 cast SequenceFile [options]
 ```

** NOTE: THE SEQUENCE FILE MUST CONTAIN PEPTIDE(S) IN FASTA FORMAT *
```
-help    ... print this text
-thr t   ... set the threshold score for reported regions
             default is 40
             t should be an integer number
-stat    ... outputs statistics information to file cast.stat
-matrix  ... use different mutation matrix (.mat) file
-verbose ... verbose mode prints filtering information to standard output
-stderr  ... verbose mode prints filtering information to standard error
```
