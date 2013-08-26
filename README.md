# ASAP2
## Automated Simultaneous Analysis Phylogenetics - Version 2

### What is it?
ASAP2 is a phylogenomic software package for identifying the evolutionary history of a set of genes. Given a set of human genes, ASAP2 will identify organisms that contain homologies for all of the human genes specified. It will then create a phylogenetic tree for each gene family (gene partition), and a consensus tree using the data contained in all gene partitions. The final consensus tree implements Partitioned Bremer Support (PBS) to predict the likelihood that the tree is, in fact, the optimal tree based on all available data.

ASAP is written primarily in the Ruby programming language, and utilizes a set of publicly-available tools to perform various data manipulations. Data input comes primarily in the form of GenBank GI numbers, with the option of including PSI-BLAST result files to further increase the set of genes to be considered in the phylogenetic analyses.

ASAP2 is the successor to ASAP - a similar tool developed by Neil Sarkar ([doi:10.1186/1471-2105-9-103](http://www.ncbi.nlm.nih.gov/pubmed/18282301)).

### System Prerequisites
Operating system: Macintosh OS X or Linux

ASAP2 requires the following software to be installed on the host machine:
* [Ruby](http://www.ruby-lang.org/en/)
    + [Bioruby](http://bioruby.org/)
    + [Terminal-table](https://github.com/visionmedia/terminal-table)
* [BLAST](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
* [TNT](http://www.zmuc.dk/public/phylogeny/tnt/)

[MUSCLE](http://www.drive5.com/muscle/) is included in the software package, and does not require separate installation by the user.

### Installation Instructions
ASAP2 can be cloned to any desired directory on a host computer. TNT requires acceptance of a license agreement, so you will have to run TNT at least once before ASAP2 has the ability to utilize it.

Before running ASAP2, you must edit the file "bin/asap2.rb" to include your email address (required for command-line queries of NCBI's GenBank and GenPept databases), and the installation directory of TNT on your computer.

### Data Workflow

* User specifies GenBank GIs of human (nucleotide) gene sequences to be analyzed
* User is given the option to include PSI-BLAST result files, prepared prior to running ASAP2
* ASAP2 performs a recursive BLASTx (nucleotide to amino-acid) on the human GIs specified previously
* The results of the recursive BLAST are combined with the PSI-BLAST result files to create a list of GenPept GIs that are highly similar to the translated human GenBank GIs being analyzed
* All of these results are used to create data partitions (see below) to be used in the phylogenetic analyses
* For each data partition, all GIs within that partition are compared against one another via BLAST-2, and the resultant expect values (E-values) are saved as easily readable charts - these charts are intended to confirm that all protein sequences are highly similar to each of the other sequences in that partition (see below regarding how to interpret these charts)
* Each data partition is aligned using MUSCLE
* Each aligned data partition is used as the basis for a phylogenetic tree specific to that partition - the resultant trees are output in parenthetical notation
* All partitions are combined into an interleaved data matrix, and the interleaved data matrix is used to create a consensus tree for all partitions, with each individual partition supplying a numerical PBS value to confirm the degrees of support for the consensus tree.

### Data File Structures

#### Results of recursive BLAST and PSI-BLAST

The GenPept sequences returned by recursive BLAST (and PSI-BLAST, if used) are stored in a file as a JSON object. The file is "results/logs/reference_sequences.json" The overall structure is as follows:

{
  "partition_1": {
    "Species_1": "gi_1",
    "Species_2": "gi_2",
    "Species_3": "gi_3"
  },
  "partition_2": {
    "Species_1": "gi_4",
    "Species_2": "gi_5",
    "Species_3": "gi_6"
  },
  "partition_3": {
    "Species_1": "gi_7",
    "Species_2": "gi_8",
    "Species_3": "gi_9"
  }
}

A sample file can be found [here](MAKE_LINK). Note that the partition names are the GI numbers of the Human GenBank sequence that was used to construct them.

#### TNT output files

INCOMPLETE

