# clustt_utils--utilities to generate sequence similarity graphs for processing and visualization

_note_: [BioLayout](http://biolayout.org) is a visualization platform for graphs; the pairs-list file can also be imported into other systems such as Cytoscape etc.; [MCL](https://www.micans.org/mcl/) is a graph clustering algorithm based on graph flow simulation, first used for sequence similarity graphs and automated protein family detection.

BioLayout:

original authors: Goldovsky L, Cases I, Enright AJ, Ouzounis CA • PMID: 16000016 DOI: 10.2165/00822942-200504010-00009

original authors: Enright AJ, Ouzounis CA • PMID: 11590107 DOI: 10.1093/bioinformatics/17.9.853

Tribe-MCL:

original authors: Enright AJ, Van Dongen S, Ouzounis CA • PMID: 11917018 DOI: 10.1093/nar/30.7.1575

original authors: Enright AJ, Kunin V, Ouzounis CA • PMID: 12888524 DOI: 10.1093/nar/gkg495

utilities by: Dimitris Vasiliou, Christos A Ouzounis [BCPLab](http://genome.academy) 2023

### mcl_clustering.sh

mcl_clustering.sh is a shell script that executes the following pipeline:

1) It accepts as input a pairs file in tabular form with 6 tabs. This input is the output of a BLAST run and has the following format:
   
   ```
   NP_085128.2-N__coords_1--1579	Hsap-XXX-01-054234	92.400	1579	120	0	1	1579	1	1579	0.02957
   NP_085128.2-N__coords_1--1579	Hsap-XXX-01-085815	92.400	1579	120	0	1	1579	1	1579	0.02956
   ```
   The script then counts the unique protein names in the first column.
   
2) It creates the [MCL](https://www.micans.org/mcl/) graph input file from the above pairs file
   by filtering only the pairs whose E-value (last column) is lower than 0.005 and
   prints the first three columns sorted by the third column (similarity score) in descending order.
   
4) Then it executes the [MCL](https://www.micans.org/mcl/) executable whose absolute path is provided as argument and uses the file generated in the previous step as input.
 
5) Finally it creates a families output file using the [MCL](https://www.micans.org/mcl/) output by prepending to each line of the [MCL](https://www.micans.org/mcl/)
   output a NR_NF string where NR is the current line (or record) number and NF the number of graph nodes in this line.
   
The script needs [MCL - a cluster algorithm for graphs](https://www.micans.org/mcl/) to be installed on the system and its executable 
path is provided to the script as an input parameter.

#### Usage:
```
Usage: mcl_clustering.sh arg1 arg2 arg3 arg4

Pairs to mcl families pipeline script

Available options:

arg1            Input tabular pairs file
arg2            Name of the MCL families output file
arg3            MCL baseline clustering parameter -I (https://www.micans.org/mcl/man/mcl.html#opt-I)
arg4            Path of the MCL executable (https://www.micans.org/mcl/)
-h, --help      Print this help and exit
-v, --verbose   Print script debug info
```
```
 ./mcl_clustering.sh input.pairs families.txt 1.8 /path_to_mcl_executable/mcl
 ```
0) The original MCL input file -> input.pairs

This creates three new files:

1) Intermediate MCL input -> input.graph
 
   This pairs-list can be used to visualize sequence similarity graphs with the [biolayout.org](http://biolayout.org) app
3) the MCL output -> out.input.graph.I1.8
4) the final families output file -> families.txt


   
