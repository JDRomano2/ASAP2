# ASAP2
## Automated Simultaneous Analysis Phylogenetics - Version 2

- - -

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

- - -

### Running ASAP2

The main script in ASAP2 is `bin/asap2.rb`. After meeting all the prereqs and configuring the TNT installation directory and default email address (see above), the script can be run by navigating to the `bin/` directory and executing:

```
ruby asap2.rb [options]
```
Under the normal (i.e., no options specified) execution procedure, ASAP2 will ask you to input the NCBI GI numbers for all of the human genes of interest, separated by spaces. For example:

```
Please enter all (nucleotide) GIs associated with the disease of interest, separated by spaces:
62088808 295844837 50083291 324021738 
```
It will then loop through each of the GIs and prompt you to either enter a PSI-BLAST output file, or leave the field blank if you are not using predetermined PSI-BLAST results for that GI. For example: 

```
Enter the relative or absolute path to PSI-BLAST results for gi 62088808 (leave blank for no PSI-BLAST file):
../data/psiblast/62088808.txt
```
ASAP2 will now perform the automated analysis. The recursive BLAST algorithm will return a message in the form of a ruby hash every time a new match is added to a data partition. It will also inform the user of how many items in the partition have been run against BLAST, along with the current total number of items in the partition. Since the data partition will likely grow as the analysis is performed, the number of total items should increment as new matches are found.

Other messages will be output to the console following the BLAST analysis. Please refer to the "Data Workflow" section below for further information on the components.

#### Options

Currently, the only option that can be specified is `resume`. In context, the command syntax is `ruby asap2.rb resume`. This option is used if you have already generated a compatible set of partitions, saved in `results/logs/blast_results.json`. This file can simply be left from a previous run of ASAP2, but is most useful when the user wants to manually create/edit the data partitions to be used for phylogenetic analysis. See below - the structure of these JSON files is enumerated in "Data file Structures - Results of recursive BLAST and PSI-BLAST".

- - -

### Data Workflow

1. User specifies GenBank GIs of human (nucleotide) gene sequences to be analyzed
2. User is given the option to include PSI-BLAST result files, prepared prior to running ASAP2
3. ASAP2 performs a recursive BLASTx (nucleotide to amino-acid) on the human GIs specified previously
4. The results of the recursive BLAST are combined with the PSI-BLAST result files to create a list of GenPept GIs that are highly similar to the translated human GenBank GIs being analyzed
5. All of these results are used to create data partitions (see below) to be used in the phylogenetic analyses
6. For each data partition, all GIs within that partition are compared against one another via BLAST-2, and the resultant expect values (E-values) are saved as easily readable charts - these charts are intended to confirm that all protein sequences are highly similar to each of the other sequences in that partition (see below regarding how to interpret these charts)
7. Each data partition is aligned using MUSCLE
8. Each aligned data partition is used as the basis for a phylogenetic tree specific to that partition - the resultant trees are output in parenthetical notation
9. All partitions are combined into an interleaved data matrix, and the interleaved data matrix is used to create a consensus tree for all partitions, with each individual partition supplying a numerical PBS value to confirm the degrees of support for the consensus tree.

- - -

### Data File Structures

#### Results of recursive BLAST and PSI-BLAST

The GenPept sequences returned by recursive BLAST (and PSI-BLAST, if used) are stored in a file as a JSON object. The file is `results/logs/reference_sequences.json` (NOTE: this file is overwritten every time you run the BLAST analysis - see "Options" for information on running the alignment and TNT analysis on a manually edited `reference_sequences.json` generated previously). The overall structure is as follows:

```
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
```

A sample file can be found [here](MAKE_LINK). Note that the partition names are the GI numbers of the Human GenBank sequence that was used to construct them.

#### TNT output files

The files generated by TNT have a couple different file formats, based on the data they contain. All of the files generated by TNT are in the `results/tnt_output` directory, and are given names based on the set of data they contain (if the data file used to create the TNT analysis is from a single gene partition, the filename will be the name of that partition followed by the appropriate filetype).

* .tre files - these are the phylogenetic trees, created in parenthetical notation. As of the current version, these can only be used in their, and must be manually converted into Newick format if you want it to be viewable in most phylogenetic tree display software - the next major release of the software will fix this, and convert the file automatically to Newick.

* .out files - These are standard TNT output files - basically a text dump of all the output generated while running the TNT script. These usually show a graphical (ASCII) representation of the resultant phylogenetic tree at the very end of the file. Note that, due to the way TNT formats output, there are often very long stretches of blank lines preceeding the trees at the end of the file, and you'll have to scroll all the way to the bottom to see the results. This may be avoidable by turning word-wrap off in the text editor you are using to view the log.