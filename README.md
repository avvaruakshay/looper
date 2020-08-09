# 100PER

[![Build](https://img.shields.io/badge/Build-passing-brightgreen.svg)]()
[![License](https://img.shields.io/badge/Licence-MIT-blue.svg)]()

<img src="./looper_logo.jpg" alt="logo" height="80"/>

100PER pronounced as Looper is a DNA tandem repeat identification tool. Tandem 
repeats are important genomic sequences which have functional and evolutionary 
significance.

Looper is scripted in C++. 

### Usage

The help message and available options can be accessed using
```bash
$ looper -h 
```
which gives the following output
```
usage: looper -i <file> -o <file> -m <int> -M <int> -l <int>

Required arguments: 
-i	<file>	Input fasta file

Optional arguments: 
-m	<int>	Minimum motif size. Default: 1
-M	<int>	Maximum motif size. Default: 6
-l	<int>	Cutoff repeat length. Default: 2*M.
 	 	Should at least be twice of maximum motif size.
-o	<file>	Output file name.Default: Input file name + _looper.tsv

```

### `-i`
**Expects:** *STRING (to be used as filename)*<br>
**Default:** *None*<br>
This is the only required argument for the program. The input file must be a 
valid FASTA file. 

### `-o`
**Expects:** *STRING (to be used as filename)*<br>
**Default:** *Input Filename + _looper.tsv (see below)*<br>
If this option is not provided, the default output filename will be the same as the input filename, with its extension replaced with '_looper.tsv'. For example, if the input filename is `my_seq.fa`, the default output filename will be `my_seq.fa_looper.tsv`. If the input filename does not have any extension, `_looper.tsv` will be appended to the filename. Please note that even in the case of no identified SSRs, the output file is still created (therefore overwriting any previous file of the same name) but with no content in the file.
#### Output for fasta
The output is a tab-delimited file, with one SSR record per line. 
The output columns follow the [BED](https://genome.ucsc.edu/FAQ/FAQformat.html) format. The details of the columns are given below:

| S.No | Column | Description |
|:----:| ------ | ----------- |
| 1 | Chromosome | Chromosome or Sequence Name as specified by the first word in the FASTA header |
| 2 | Repeat Start | 0-based start position of SSR in the Chromosome |
| 3 | Repeat Stop | End position of SSR in the Chromosome |
| 4 | Repeat Class | Class of repeat as grouped by their cyclical variations |
| 5 | Repeat Length | Total length of identified repeat in nt |
| 6 | Repeat Strand | Strand of SSR based on their cyclical variation |
| 7 | Motif Number | Number of times the base motif is repeated |
| 8 | Actual Repeat | Starting sequence of the SSR irrespective of Repeat class and strand|

An example output showing some of the largest repeats from *Drosophila melanogaster* is given below
```
X       22012826  22014795  ACTGGG  1969    -       328     TCCCAG
2RHet   591337    591966    AATACT  629     -       104     ATTAGT
4       1042143   1042690   AAATAT  547     +       91      AAATAT
2RHet   598244    598789    AATACT  545     -       90      AGTATT
XHet    122       663       AGAT    541     +       135     GATA
X       22422335  22422827  AGAT    492     +       123     GATA
3R      975265    975710    AAAT    445     -       111     TTAT
X       15442288  15442724  ACAGAT  436     +       72      ACAGAT
2L      22086818  22087152  AATACT  334     -       55      TATTAG
YHet    137144    137466    AAGAC   322     -       64      CTTGT
```

### `-m`
**Expects:** *INTEGER*<br>
**Default:** *1*<br>
Minimum length of motifs to be considered. By default, looper ignores redundant 
motifs. For example, a stretch of 12 A's is considered a monomer repeat of 12 
A's rather than a dimer repeat of 6 AA's. 

### `-M`
**Expects:** *INTEGER*<br>
**Default:** *6*<br>
Maximum length of motifs to be considered. Setting a large value of `-M` has a 
non-trivial effect on both the runtime and memory usage of looper.

### `-l`
**Expects:** *INTEGER*<br>
**Default:** *2* * *M*<br>
Minimum length cut-off to be considered when finding an SSR. The same cut-off 
will apply for SSRs of all motif lengths, even if the motif length is not a 
divisor of this value. In such cases, SSRs that end with a partial motif are 
also picked if they pass the length cut-off. This value should be at least twice
of the maximum motif size.

## Contact
For queries or suggestions, please contact:

Divya Tej Sowpati - <tej@ccmb.res.in><br>
Akshay Kumar Avvaru - <avvaru@ccmb.res.in>