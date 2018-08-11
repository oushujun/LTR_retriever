### Introduction ###

LTR_retriever is a command line program (in Perl) for accurate identification of LTR retrotransposons (LTR-RTs) from outputs of LTRharvest, LTR_FINDER, and/or MGEScan-LTR and generating non-redundant LTR-RT library for genome annotations.

By default, the program will generate whole-genome LTR-RT annotation and the LTR Assembly Index (LAI) for evaluations of the assembly continuity of the input genome. Users can also run LAI separately (see `Usage`).

### Installation ###

To run LTR_retriever you need to provide the paths to the following dependent programs.
1. makeblastdb, blastn, and blastx in the BLAST+ package (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/),
2. cd-hit-est in the CDHIT package (http://weizhongli-lab.org/cd-hit/) OR 
   blastclust in the BLAST package (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy/2.2.25/),
3. hmmsearch in the HMMER package (http://hmmer.org/), and
4. RepeatMasker (http://www.repeatmasker.org/).

Simply modify the 'paths' file in the LTR_retriever directory

	vi /your_path_to/LTR_retriever/paths

Then modify lines below:

	BLAST+=/your_path_to/BLAST+2.2.30/bin/
	RepeatMasker=/your_path_to/RepeatMasker4.0.0/
	HMMER=/your_path_to/HMMER3.1b2/bin/
	CDHIT=/your_path_to/CDHIT4.6.1/
	BLAST=/your_path_to/BLAST2.2.26/bin/ #not required if CDHIT provided


### Inputs ###

Two types of inputs are needed for LTR_retriever
1. Genomic sequence
2. LTR-RT candidates

LTR_retriever takes multiple LTR-RT candidate inputs including the screen output of LTRharvest, the screen output of LTR_FINDER, and the candidate output of MGEScan-LTR. Users need to obtain the input file(s) from the aforementioned programs before running LTR_retriever. Either a single input source or a combination of multiple inputs are acceptable. For more details and examples please see the manual.

### Outputs ###

The output of LTR_retriever includes:
1. Intact LTR-RTs with coordinate and structural information
	- Summary tables (.pass.list)
	- GFF3 format output (.pass.list.gff3)
2. LTR-RT library
	- All non-redundant LTR-RTs (.LTRlib.fa)
	- All non-TGCA LTR-RTs (.nmtf.LTRlib.fa)
	- All LTR-RTs with redundancy (.LTRlib.redundant.fa)
3. Whole-genome LTR-RT annotation by the non-redundant library
	- GFF format output (.out.gff)
	- LTR family summary (.out.fam.size.list)
	- LTR superfamily summary (.out.superfam.size.list)
	- LTR distribution on each chromosome (.out.LTR.distribution.txt)
4. LTR Assembly Index (.out.LAI)

### Usage ###

To run LTR_retriever:

	/your_path_to/LTR_retriever -genome genomefile -inharvest LTRharvest_input [options]

To run LAI:

	/your_path_to/LAI -genome genome.fa -intact intact.pass.list -all genome.out [options]

For more details about the usage and parameter settings, please see the help pages by running:

	/your_path_to/LTR_retriever -h

	/your_path_to/LAI -h
	
Or refer to the manual document.


For questions and Issues Please See: https://github.com/oushujun/LTR_retriever/issues

### Citations ###

For LTR_retriever, please cite:

`Ou S. and Jiang N. (2018). LTR_retriever: A Highly Accurate and Sensitive Program for Identification of Long Terminal Repeat Retrotransposons. Plant Physiol. 176(2): 1410-1422.`

For LAI, please cite:

`Ou S., Chen J. and Jiang N. (2018). Assessing genome assembly quality using the LTR Assembly Index (LAI). Nucleic Acids Res. gky730:` https://doi.org/10.1093/nar/gky730
