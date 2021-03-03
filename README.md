# Looper

[![Build](https://img.shields.io/badge/Build-passing-brightgreen.svg)]()
[![License](https://img.shields.io/badge/Licence-MIT-blue.svg)]()

<img src="./lib/looper_logo.jpg" alt="logo" height="80"/>

Looper is a DNA tandem repeat identification tool. Tandem repeats are important 
genomic sequences which have functional and evolutionary significance.

Looper is scripted in C++. 

### Usage

The help message and available options can be accessed using
```bash
$ python looper.py -h 
```
which gives the following output
```
usage: looper.py [-h] -i <FILE> [-o <FILE>] [-m <INT>] [-M <INT>] [-l <INT>] [-a] [-g <FILE>]
                 [--anno-format <STR>] [--gene-key <STR>] [--up-promoter <INT>]
                 [--down-promoter <INT>]

Required arguments:
  -i, --input	<FILE>	Input sequence file. Fasta format.

Optional arguments:
  -o, --output			<FILE>	Output file name. Default: Input file name + _looper.tsv
  -m, --min-motif-size	<INT>	Minimum size of a repeat motif in bp. Default: 1
  -M, --max-motif-size	<INT>	Maximum size of a repeat motif in bp. Default: 6
  -l, --min-length		<INT>	Cutoff repeat length. Default: 2*M.
 								Should at least be twice of maximum motif size.
  -a, --analyse					Generate a summary HTML report.

Compound repeat arguments:
  --compound            Report compound repeats. The output of compound repeats is a separate file with the suffix ".compound".
  -d <INT>, --comp-dist <INT>
                        Maximum distance between individual repeats of compound repeat. Use negative to denote overlap. Default: 0

Annotation arguments:
  -g, --annotate		<FILE>	Genomic feature file to annotate repeats w.r.t genes.
  								Both GFF and GTF can be processed.
  --anno-format			<STR>   Format of genomic feature file. 
  								Valid inputs: GFF, GTF. Default: GFF
  --gene-key			<STR>   Attribute used as unique identifier for gene name.
  								The default identifier is "gene". 
  --up-promoter			<INT>   Upstream distance(bp) from TSS to be considered as
  								promoter region. Default: 1000
  --down-promoter		<INT>   Downstream distance(bp) from TSS to be considered as
  								promoter region. Default: 1000
```

### `-i or --input`
**Expects:** *STRING (to be used as filename)*<br>
**Default:** *None*<br>
This is the only required argument for the program. The input file must be a 
valid FASTA file. 

### `-o or --output`
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

### `-m or --min-motif-size`
**Expects:** *INTEGER*<br>
**Default:** *1*<br>
Minimum length of motifs to be considered. By default, looper ignores redundant 
motifs. For example, a stretch of 12 A's is considered a monomer repeat of 12 
A's rather than a dimer repeat of 6 AA's. 

### `-M or --max-motif-size`
**Expects:** *INTEGER*<br>
**Default:** *6*<br>
Maximum length of motifs to be considered. Setting a large value of `-M` has a 
non-trivial effect on both the runtime and memory usage of looper.

### `-l or --min-length`
**Expects:** *INTEGER*<br>
**Default:** *2* * *M*<br>
Minimum length cut-off to be considered when finding an SSR. The same cut-off 
will apply for SSRs of all motif lengths, even if the motif length is not a 
divisor of this value. In such cases, SSRs that end with a partial motif are 
also picked if they pass the length cut-off. This value should be at least twice
of the maximum motif size.

### `-a or --analyze`
**Expects:** *None*<br>
**Default:** *False*<br>
In addition to the default tab-separated output, Looper can also generate a fully
interactive HTML report for easy downstream analysis of the repeat data. The 
filename will be the same prefix as that of the main output. For example, if the
input filename was my_seq.fa, the analysis report will be my_seq_looper.html. An 
example HTML report, generated from the repeat data of Homo sapiens (build hg19),
can be accessed here (Right click -> Save As).

### `--compound`
**Expects:** *None*<br>
**Default:** *False*<br>
This is flag which when set to true reports all compound repeats. Compound repeats
are repeats which are either overlapping or separated by a small gap. Compound repeats
are reported in a separate file with the extension ".compound". The maximum distance
between two individual repeats of a compound repeat can be set using the option `-d`.
The compoud repeat has four columns with first four denoting the sequence, start, 
end and length of the repeat. The fifth column denotes the repeat classes of individual
repeats of the compound repeat and the number of different cyclical variations that
have occured as compound repeat. The sixth column is the strand of the individual
repeats. The seventh column is denotes the actual motifs of individual repeats, 
repeat length and the distance between the individual repeats. Below is an example
output reporting compound repeats.

```
chr1	10330	10392	(AACCCT)2	(CCCTAA)23-(0)-(ACCCTA)40	+|+
chr1	10397	10468	(AACCCT)2	(CCCTAA)45-(0)-(CCCTAA)28	+|+
chr1	10630	10654	(CCGCG)1(AGGCGC)1	(CGCGC)12-(0)-(GGCGCA)13	+|+
chr1	10659	10683	(CCGCG)1(AGGCGC)1	(CGCGC)12-(0)-(GGCGCA)13	+|+
chr1	10688	10712	(CCGCG)1(AGGCGC)1	(CGCGC)12-(0)-(GGCGCA)13	+|+
chr1	10717	10741	(CCGCG)1(AGGCGC)1	(CGCGC)12-(0)-(GGCGCA)13	+|+
chr1	10746	10770	(CCGCG)1(AGGCGC)1	(CGCGC)12-(0)-(GGCGCA)13	+|+
chr1	10775	10799	(CCGCG)1(AGGCGC)1	(CGCGC)12-(0)-(GGCGCA)13	+|+
chr1	10847	10871	(CCGCG)1(AGGCGC)1	(CGCGC)12-(0)-(GGCGCA)13	+|+
chr1	27543	27561	(AAAAG)1(AAAAAG)1	(TTTCT)13-(0)-(TTTTCT)14	-|-
chr1	49834	49861	(AAAC)1(AAAAC)1	(AAAC)19-(0)-(AAACA)15	+|+
chr1	50560	50587	(AAAAG)1(AAAGG)1	(TTCTT)12-(0)-(TTTCC)17	-|-
chr1	66257	66287	(AAATAT)1(AATAT)2	(AATATA)13-(0)-(AATAT)12-(-2)-(ATTAT)12	+|+|-
chr1	66294	66312	(AAATAT)1(AATAT)1	(AATATA)13-(0)-(AATAT)12	+|+
```

### `-d or --comp-dist`
**Expects:** *INTEGER*<br>
**Default:** *0*<br>
Maximum distance between two individual repeats of a compound repeat. To report repeats overlapping by X bp use negative numbers i.e., -X. The default is 0bp. 

### `-g or --annotate`
**Expects:** *FILE*<br>
**Default:** *None*<br>
Input a genomic feature file to annotate the repeats in the genomic context. 
Looper accepts both GFF and GTF format genomic feature files. Each repeat is 
annotated w.r.t the closest gene and classified either as Genic, Exonic, 
Intronic and Intergenic according to the position of the repeat. Besides this, 
the repeat is also checked if it falls in the promoter region of the gene. 
Annotation adds 7 columns to the default looper output which already consist 8 
columns.

| S.No | Column | Description |
|:----:| ------ | ----------- |
| 9 | Gene name | Name of the closest gene |
| 10 | Gene Start | Start position of gene in the Chromosome |
| 11 | Gene Stop | End position of gene in the Chromosome |
| 12 | Strand | The strand orientation of the gene |
| 13 | Genomic annotation | Annotation of the repeat w.r.t to the gene. Possible annotations are {Genic, Exonic, Intronic, Intergenic} |
| 14 | Promoter annotation | If repeat falls in the promoter region of the closest gene. The default promoter region is 1Kb upstream and downstream of TSS. |
| 15 | Distance from TSS | Distance of the repeat from the TSS of the gene. |

### `--anno-format`
**Expects:** *STRING*<br>
**Default:** *GFF*<br>
Option to specify the format of the input genomic feature file. Accepted inputs 
are GFF or GTF. More details about the GFF and GTF formats can be found 
[here](https://asia.ensembl.org/info/website/upload/gff.html).

### `--gene-key`
**Expects:** *STRING*<br>
**Default:** *gene*<br>
The attribute key used for the name of the gene in the GFF/GTF file. In the 
below example GFF file, we have the location of a gene and it's mRNA and exon 
locations. The last column of the file specifies attributes associated with each
feature, like ID, Parent, gene etc. Looper uses on of the attribute to identify 
the gene and also it's exons. In th below example the key "gene" can be used to 
identify gene and the exons of the gene as they have the same gene name. Please 
check your GFF/GTF file for a robust attribute key which can identify all genes 
and their corresponding exons. We are actively working on better annotation 
where we can identify genes and their exons based on the ID and Parent.

```
# Sample GFF
NC_004354.4	RefSeq	gene	124370	126714	.	-	.	ID=gene1;Name=CG17636;gbkey=Gene;gene=CG17636;gene_biotype=protein_coding;gene_synonym=DmelCG17636,EG:23E12.1;
NC_004354.4	RefSeq	mRNA	124370	126714	.	-	.	ID=rna1;Parent=gene1;Name=NM_001103384.3;gbkey=mRNA;gene=CG17636;transcript_id=NM_001103384.3
NC_004354.4	RefSeq	exon	126626	126714	.	-	.	ID=id13;Parent=rna1;gbkey=mRNA;gene=CG17636;transcript_id=NM_001103384.3
NC_004354.4	RefSeq	exon	125495	126259	.	-	.	ID=id14;Parent=rna1;gbkey=mRNA;gene=CG17636;transcript_id=NM_001103384.3
```

### `--up-promoter`
**Expects:** *INT*<br>
**Default:** *1000*<br>
Upstream distance(bp) from the TSS of the gene to be considered as promoter 
region. Default 1000.

### `--down-promoter`
**Expects:** *INT*<br>
**Default:** *1000*<br>
Downstream distance(bp) from the TSS of the gene to be considered as promoter 
region. Default 1000.

## Contact
For queries or suggestions, please contact:

Divya Tej Sowpati - <tej@ccmb.res.in><br>
Akshay Kumar Avvaru - <avvaru@ccmb.res.in>