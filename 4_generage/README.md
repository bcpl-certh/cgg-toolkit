# generage--a robust algorithm for sequence clustering and domain detection

Generage:

original authors: Enright AJ, Ouzounis CA • PMID: 10871267 DOI: 10.1093/bioinformatics/16.5.451 

updates by: Dimitris Vasiliou, Christos A Ouzounis, [BCPLab](http://genome.academy/) 2023

DifFuse:

original authors: Enright AJ, Iliopoulos I, Kyrpides NC, Ouzounis CA • PMID: 10573422 DOI: 10.1038/47056

updates by: Dimitris Vasiliou, Christos A Ouzounis, [BCPLab](http://genome.academy/) 2023

#### Generage usage:
```
export PRSSDIR=$(pwd)

make

./generage file1 file2
```
where file1 is a FASTA format protein database and file2 is a parsed BLAST results file.

#### DifFuse usage:
```
cd diffuse

make

./genefuse meth.fasta sach2.fasta meth_v_meth.m8.blastout meth_v_sach.m8.blastout -verbose -rand 1000
```

** NOTE: See generage.sh and diffuse/genefuse.sh for actual call examples of the generage and genefuse executables **

#### Dependencies:

#### [fasta-35.4.12](https://fasta.bioch.virginia.edu/wrpearson/fasta/)

Copyright 1988, 1991, 1992, 1993, 1994 1995, by William
R. Pearson and the University of Virginia.  All rights
reserved. The FASTA program and documentation may not be sold or
incorporated into a commercial product, in whole or in part,
without written consent of William R. Pearson and the University
of Virginia.  For further information regarding permission for
use or reproduction, please contact:

David Hudson
Assistant Provost for Research
University of Virginia
P.O. Box 400301
Charlottesville, VA  22906-9025

(434) 924-3606
     
Code in the smith_waterman_sse2.c and smith_waterman_sse2.h files
is copyright (c) 2006 by Michael Farrar.

This program may not be sold or incorporated into a commercial
product, in whole or in part, without written consent of Michael
Farrar.  For further information regarding permission for use or
reproduction, please contact: Michael Farrar at
farrar.michael@gmail.com.


