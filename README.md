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