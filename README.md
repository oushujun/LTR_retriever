### Introduction ###

LTR_retriever is a command line program (in Perl) for accurate identification of LTR retrotransposons (LTR-RTs) from outputs of LTRharvest, LTR_FINDER, and/or MGEScan-LTR and generating non-redundant LTR-RT library for genome annotations.

### Installation ###

To run LTR_retriever you need to provide the paths to the following dependent programs.
1. makeblastdb, blastn, and blastx in the BLAST+ package (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/),
2. cd-hit-est in the CDHIT package (http://weizhongli-lab.org/cd-hit/) OR 
   blastclust in the BLAST package (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy/2.2.25/),
3. hmmsearch in the HMMER package (http://hmmer.org/), and
4. RepeatMasker (http://www.repeatmasker.org/).

Simply modify the 'paths' file in the same directory
In your terminal, type:

	vi /your_path_to/LTR_retriever/paths

Then modify lines below:

	BLAST+=/your_path_to/BLAST+2.2.30/bin/
	RepeatMasker=/your_path_to/RepeatMasker4.0.0/
	HMMER=/your_path_to/HMMER3.1b2/bin/
	CDHIT=/your_path_to/CDHIT4.6.1/
	BLAST=/your_path_to/BLAST2.2.26/bin/ #not required if CDHIT provided`


### Inputs ###

Two types of inputs are needed for LTR_retriever
1. Genomic sequence
2. LTR-RT candidates

LTR_retriever takes multiple LTR-RT candidate inputs including the screen output of LTRharvest, the screen output of LTR_FINDER, and the candidate output of MGEScan-LTR. Users need to obtain the input file(s) from the aforementioned programs before running LTR_retriever. Either a single input source or a combination of multiple inputs are acceptable. For more details and examples please see the manual.

### Outputs ###

The output of LTR_retriever includes:
1. A summary table for the identified intact LTR-RTs with coordinate and structural information (.pass.list)
2. A non-redundant LTR-RT library (exemplar) in the FASTA format (.LTRlib.fa)
3. A redundant LTR-RT library in the FASTA format (.LTRlib.redundant.fa)
4. A GFF3 format file for all intact LTR-RTs (.pass.list.gff3)
5. A GFF format file for the whole-genome LTR-RT annotation (.gff)
6. A table for the LTR Assembly Index (.LAI)

### Usage ###

	LTR_retriever -genome genomefile -inharvest LTRharvest_input [options]

For more details about the usage and parameter settings, please see the help page of LTR_retriever by running:
	./LTR_retriever -h
Or refer to the manual document.

For questions and Issues Please See: https://github.com/oushujun/LTR_retriever/issues
