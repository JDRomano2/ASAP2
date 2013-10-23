#!/usr/bin/env ruby

#NEXT: FIX INVALID STRINGS IN NUC_REFERENCE_SEQUENCES.JSON!

################################################################################
#
# asap2.rb
#
# Copyright: Copyright (C) 2013 by Joe Romano
# Contact: jdromano2@gmail.com
# License: GNU General Public License
#
################################################################################
#
# File Description:
# => ASAP2 - Automated Simultaneous Analysis Phylogenetics; Ver. 2
# => Still in alpha development phase
#
################################################################################

require 'rubygems'
#require 'bundler/setup'
#Bundler.setup
#Bundler.require
require 'json'
require 'pp'
require 'bio'
require 'terminal-table'

################################################################################
######################
# USER CONFIGURATION #
######################

$user_email = "jdromano@uvm.edu"

$tnt_path = "/usr/local/bin/tnt"

##########################
# END USER CONFIGURATION #
##########################
################################################################################

puts "Run full program or resume with TNT analysis of BLAST results?\nFull program - 1\nResume analysis - 2\n"
opt1 = gets.chomp()
if opt1 == "1"
  resume = 0
elsif opt1 == "2"
  resume = 1
else
  puts "Invalid selection - proceeding with full analysis."
  resume = 0
end

=begin
puts "Run analysis on nucleotide or protein results?\nNucleotide - 1\nProtein - 2\n"
opt2 = gets.chomp()
if opt2 == "1"
  mode = 1
elsif opt2 == "2"
  mode = 2
else
  puts "Invalid selection - proceeding in Protein mode."
  mode = 2
end
=end
mode = 2 #Only allow protein mode for now - may fix nucleotide mode later...

#Set up REST and local parameters
Bio::NCBI.default_email = $user_email
ncbi_fetch = Bio::NCBI::REST::EFetch.new()

class Get_NCBI_Sequence
  def initialize(accessor, sequence = '')
    @accessor = accessor
    @sequence = sequence
  end

  def fetch_sequence
    ncbi_fetch = Bio::NCBI::REST::EFetch.new()
    @sequence = Bio::FastaFormat.new(ncbi_fetch.sequence(@accessor, "fasta"))
    return @sequence
  end

  attr_accessor :sequence
end

class Splice
  attr_accessor :fasta
  attr_accessor :accession
  attr_accessor :seq_spliced

  def initialize(complement, source)
    @ncbi_fetch = Bio::NCBI::REST::EFetch.new()
    @complement = complement
    @source = source
    @organism = @source.tr(" ", "_")
    if (@organism == "")
      raise "Passed blank organism name to Splice.initialize"
    end
  end

  def parse_components
    @accession = @complement.scan(/^.*.{11}?join\(([0-9A-Z]+)\..*$/)
    @ranges = @complement.scan(/:(\d+\.\.\d+)/)
  end

  def make_fasta
    @seq = Bio::Sequence.auto(@ncbi_fetch.sequence(@accession, "fasta").gsub!(/^>.+\n/, '')) #gsub removes header line from fasta
    @range_string = @ranges.join(",")
    if (@complement.scan(/^.*.{11}join\(([0-9A-Z]+)\..*$/).length > 0)
      @seq_spliced = @seq.splice("complement(join(#{@range_string}))")
    else
      @seq_spliced = @seq.splice("join(#{@range_string})")
    end
    @fasta = Bio::FastaFormat.new(">#{@organism} transsplice mRNA coding protein #{@accession}\n#{@seq_spliced}")
  end
end

# Set global variables
case
when RUBY_PLATFORM.downcase.include?("linux")
  $platform = :linux
when RUBY_PLATFORM.downcase.include?("darwin")
  $platform = :darwin
end
$pwd = `pwd`

# Validate dependencies: BLAST, Muscle, TNT
begin
  $blastx = `which blastx`
  # Handle exceptions? after next line too
  $blastp = `which blastp`
  $blastn = `which blastn`
  # Muscle prepackaged - doesn't require user configuration
  $platform == :linux ? $muscle = "muscle3.8.31_i86linux64" : nil
  $platform == :darwin ? $muscle = "./muscle3.8.31_i86darwin64" : nil
  $platform == :linux ? $tnt = "tnt" : nil
  $platform == :darwin ? $tnt = "./tnt.command" : nil
rescue Exception
  $stderr.print "Couldn't validate dependencies\n"
  raise
end

################################################################################
# BEGIN BLAST ANALAYSIS OF SPECIFIED SEQUENCES
################################################################################

if resume == 0 #Without "resume" option, proceed with recursive BLAST + PSI-BLAST as normal
  # Collect GIs of target disease; store in an array
  puts "Please enter all (nucleotide) GIs associated with the disease of interest, separated by spaces:"
  human_gis = gets.split(" ")

  # Make hash that will contain results of all BLASTs performed
  # Refer to ../doc/data_structure.txt for further info
  blast_results = Hash.new()
  # Also, define the names of the keys (corresponding to each user-entered GI), to avoid problems down the road:
  human_gis.each { |name| blast_results[name] = Hash.new() }

  # Get PSI-BLAST filenames one-by-one, if user wants to specify them
  all_psiblast_filenames = Hash.new()
  human_gis.each do |gi|
    puts "Enter the filename of PSI-BLAST results for gi #{gi} in ../data/psiblast_files (leave blank for no PSI-BLAST file): "
    curr_psiblast = gets.chomp
    all_psiblast_filenames[gi] = curr_psiblast
  end

  all_organisms_represented = Array.new()

  # Iterate through array, running recursive BLAST on each, store results in a hash, each key is a gene and each index is a hash
  human_gis.each do |gis|

    puts "Running recursive blast on gi #{gis}"

    psiblast_filename = all_psiblast_filenames[gis]

    #instantiate variables
    current_gi = String.new()
    current_taxon = String.new()
    gi_list = Array.new()
    used_gis = Array.new() #Needed to prevent duplicates
    psiblast_gis = Array.new()
    both = Array.new()
    psiblast_only = Array.new()
    recursive_blast_only = Array.new()
    organisms_represented = Array.new()

    #Store query GI sequence as FastaFormat object and in temporary local file
    query_fasta = Bio::FastaFormat.new(ncbi_fetch.sequence(gis, format = 'fasta'))
    File.open('blast_input.fasta', 'w') {|f| f.puts(query_fasta)}

    #Open psiblast results file and save as string
    #Note: the .insert() method takes care of relative filepath
    psiblast_filename == "" ? psiblast_results = "" : psiblast_results = File.open(psiblast_filename.insert(0, "../data/psiblast_files/").chomp, 'r') { |f| f.read }

    #Run appropriate BLAST (depending on analysis type) on query gi and store results as string (tab-delimited)
    if mode == 2
      puts "Running BLASTx on query gi. This may take a few minutes.\n"
      query_blast = `blastx -remote -db nr -evalue 1e-256 -outfmt 6 -query blast_input.fasta`
    elsif mode == 1
      puts "Running BLASTn on query gi. This may take a few minutes.\n"
      query_blast = `blastn -remote -db nr -evalue 1e-256 -outfmt 6 -query blast_input.fasta`
    end

    #Parse the GI from each BLASTx result, and store that GI and associated attributes in hashes
    #Add each of these hashes to the array "gi_list"
    query_blast.each_line do |line|
      line.scan(/^gi.*?(gi\|\d+).*?$/) {|gi| current_gi = "#{gi.join}"}
      unless used_gis.include?(current_gi)
        current_sequence = Bio::GenPept.new(ncbi_fetch.sequence(current_gi))
        current_organism = current_sequence.organism()
        organisms_represented.push(current_organism) unless organisms_represented.include?(current_organism)
        ncbi_fetch.sequence(current_gi).scan(/\/db_xref="taxon:(\d+)"/) {|taxon| current_taxon = "#{taxon.join}"}
        used_gis.push(current_gi)
        current_hash = { :gi => current_gi, :organism => current_organism, :taxon => current_taxon }
        gi_list.push(current_hash)
        puts "Pushed #{current_hash}"
      end
    end
    puts "\n"

    #Iterate through each entry, running appropriate BLAST and appending any new records to
    #growing array of hashes gi_list
    gi_list.each_with_index do |record, index|
      puts "Running BLAST on record \##{index + 1} of #{gi_list.length}"
      puts "Current GI is #{record[:gi]}"
      current_fasta = Bio::FastaFormat.new(ncbi_fetch.sequence(record[:gi], format = 'fasta'))
      File.open('blast_input.fasta', 'w') {|f| f.puts(current_fasta)}
      if mode == 2
        current_blast = `blastp -remote -db nr -evalue 1e-256 -outfmt 6 -query blast_input.fasta`
      elsif mode == 1
        current_blast = `blastn -remote -db nr -evalue 1e-256 -outfmt 6 -query blast_input.fasta`
      end
      current_blast.each_line do |line|
        line.scan(/^gi.*?(gi\|\d+).*?$/) {|gi| current_gi = "#{gi.join}"}
        unless used_gis.include?(current_gi)
          current_sequence = Bio::GenPept.new(ncbi_fetch.sequence(current_gi))
          current_organism = current_sequence.organism()
          organisms_represented.push(current_organism) unless organisms_represented.include?(current_organism)
          ncbi_fetch.sequence(current_gi).scan(/\/db_xref="taxon:(\d+)"/) {|taxon| current_taxon = "#{taxon.join}"}
          used_gis.push(current_gi)
          current_hash = { :gi => current_gi, :organism => current_organism, :taxon => current_taxon }
          gi_list.push(current_hash)
          puts "Pushed #{current_hash}"
        end
      end
    end

    #Iterate through lines of psiblast results and pull GIs - add new GIs to an
    #array
    current_gi = String.new()
    psiblast_results.each_line do |line|
      line.scan(/^gi.*?(gi\|\d+).*?$/) {|gi| current_gi = "#{gi.join}"}
      current_sequence = Bio::GenPept.new(ncbi_fetch.sequence(current_gi))
      current_organism = current_sequence.organism()
      organisms_represented.push(current_organism) unless organisms_represented.include?(current_organism)
      ncbi_fetch.sequence(current_gi).scan(/\/db_xref="taxon:(\d+)"/) {|taxon| current_taxon = "#{taxon.join}"}
      current_hash = { :gi => current_gi, :organism => current_organism, :taxon => current_taxon }
      psiblast_gis.push(current_hash)
    end

    #Sort all GIs into arrays, specifying whether the GI is present in either
    #the PSI-BLAST results, the recursive BLAST results, or in both
    gi_list.each do |recursive_blast_record|
      if psiblast_gis.include?(recursive_blast_record)
        both.push(recursive_blast_record)
      else
        recursive_blast_only.push(recursive_blast_record)
      end
    end
    psiblast_gis.each do |psiblast_record|
      unless gi_list.include?(psiblast_record)
        psiblast_only.push(psiblast_record)
      end
    end

    # print the results for the recursive BLAST to STDOUT
    puts "\n\n\n###GIs in both sets:\n"
    pp both
    puts "###"
    puts "size: #{both.length} items"
    puts "\n\n\n###GIs only in recursive BLAST:\n"
    pp recursive_blast_only
    puts "###"
    puts "size: #{recursive_blast_only.length} items"
    puts "\n\n\n###GIs only in PSI-BLAST:\n"
    pp psiblast_only
    puts "###"
    puts "size: #{psiblast_only.length} items"

    # find intersection of all three lists
    items_flattened = both | recursive_blast_only | psiblast_only

    # Create a hash for each set of values, and store it in
    blast_results[gis] = items_flattened

    all_organisms_represented.push(organisms_represented)
  end
  File.open('../results/logs/blast_results.json', 'w') {|f| f.puts(JSON.pretty_generate(blast_results))}

  # find species present in all partitions
  model_species = all_organisms_represented.inject(:&)

  # CREATE HASH WITH REFERENCE SEQUENCE FOR EACH POTENTIAL MODEL BY PARTITION
  trimmed_results = Hash.new()
  blast_results.each do |partition, hits|
    gis_per_partition = Hash.new()
    model_species.each do |species|
      current = String.new()
      i = 0 #increment through array of hashes
      until current.length > 0 do
        if hits[i][:organism] == species
          current = hits[i][:gi]
        end
        i += 1
      end
      gis_per_partition[species] = current
    end
    trimmed_results[partition] = gis_per_partition
  end
  pp "\n"
  pp trimmed_results
  File.open('../results/logs/reference_sequences.json', 'w') {|f| f.puts(JSON.pretty_generate(trimmed_results))}
end
#NOTE: Removed conditional test for results! This homogenizes the name
#"results" instead of "trimmed_results" - revert if this creates issues...
results = JSON.parse( IO.read('../results/logs/reference_sequences.json') )

################################################################################
# FIND NUCLEOTIDE SEQUENCES ENCODING PROTEIN SEQUENCES
################################################################################

=begin
# Determine nucleotide sequences that encode the protein reference sequences
nucleotide_results = results.clone
if true
  results.each do |humangi, partition|
    partition.each do |species, gi|
      puts "Finding nucleotide sequence encoding protein #{gi}"
      record = Bio::GenPept.new(ncbi_fetch.sequence(gi))
      cds = record.features[(record.features.length - 1)]
      codedby = nil
      cds.qualifiers.each { |a_qualifier| codedby = a_qualifier.value if a_qualifier.qualifier == "coded_by" }
      codedby.sub!(/\..*$/, '') unless codedby.scan(/^.{11}?join\(.+/).length > 0
      if codedby == nil
        raise ArgumentError.new("Could not find nucleotide sequence for #{gi}. You should look the protein sequence up on www.ncbi.nlm.nih.gov, and find the reference sequence. Then, in 'reference_sequences.json', change #{gi} to the correct reference sequence, and restart this script, selecting option #2 at the first prompt.")
      elsif (codedby.scan(/^.{11}?join\(.+/).length > 0)
        puts "found transsplice gene #{codedby} - expect a long string in logfile!"
      elsif (codedby.match(/^[A-Z]{2}_?\d{4,}$/) == nil && codedby.match(/^[A-Z]_?\d{5,}$/) == nil)
        puts "codedby = '#{codedby}'"
        raise ArgumentError.new("Invalid CDS data for protein sequence #{gi}. Sequence is likely large genomic construct. Look up protein sequence and find a better alternative.")
      end
      nucleotide_results[humangi][species] = codedby
    end
  end
end
File.open('../results/logs/nuc_reference_sequences.json', 'w') {|f| f.puts(JSON.pretty_generate(nucleotide_results))}
=end

nucleotide_results = JSON.parse( IO.read('../results/logs/nuc_reference_sequences.json'))

################################################################################
# CREATE TERMINAL TABLES FOR IDENTIFIED PROTEIN SEQUENCES
################################################################################
=begin
all_tables = Array.new()
all_tables.push("(Remember to turn word-wrap off, or tables will not be readable on most displays)")
# create the first row of each partition's terminal-table - this is static for all partitions
firstrow = Array.new()
firstrow.push("") # row 1, col 1 is blank
results.first[1].each_key { |organism| firstrow.push(organism) }
# loop through each partition, creating a terminal table for each
results.each do |gene, matches|
  puts "Expect value table for #{gene}"
  rows = Array.new()
  rows.push(firstrow)
  firstrow.drop(1).each do |organism|
  puts "organism is #{organism}"#
    #build rows 2-...
    curr_row = Array.new()
    curr_row.push(organism)
    subject = results[gene][organism]
    subject_fasta = Bio::FastaFormat.new(ncbi_fetch.sequence(subject, format = 'fasta'))
    File.open('blast_subject.fasta', 'w') {|f| f.puts(subject_fasta)}
    firstrow.drop(1).each do |organism2|
      query = results[gene][organism2]
      query_fasta = Bio::FastaFormat.new(ncbi_fetch.sequence(query, format = 'fasta'))
      File.open('blast_query.fasta', 'w') {|f| f.puts(query_fasta)}
      if mode == 2
        current_blast = `blastp -query blast_query.fasta -subject blast_subject.fasta -outfmt '6 evalue'`
      elsif mode == 1
        current_blast = `blastn -query blast_query.fasta -subject blast_subject.fasta -outfmt '6 evalue'`
      end
      evalue = current_blast.lines.first
      curr_row.push(evalue)
    end
    rows.push(curr_row)
  end
  # Terminal Table-specific syntax - refer to Terminal Table documentation
  table = Terminal::Table.new :rows => rows
  puts table # Will display terminal tables to user before writing to file
  puts "\n\n"
  all_tables.push(table) # Store all terminal tables in array
end
#Can't print terminal-table objects to string by default
#Instead, temporarily redirect STDOUT to file, and print the terminal-tables
orig_std_out = STDOUT.clone #make clone of STDOUT to restore normal functionality later
STDOUT.reopen(File.open('../results/logs/eval_tables.txt', 'w')) #instruct where to pipe STDOUT
all_tables.each { |table| puts table; puts "\n" } #print each table
STDOUT.reopen(orig_std_out) #switch back to normal STDOUT
=end

################################################################################
# FIRST, MAKE TREES FOR PROTEIN SEQUENCES...
################################################################################

=begin
################################################################################
# CREATE AND ALIGN FASTA FILES
################################################################################
# Create a fasta file for each partition, pulling sequences for all GIs from NCBI
results.each_pair do |partition, sequences|
  fasta_contents = String.new()
  sequences.each do |gi|
    pp gi
    fasta_contents << Get_NCBI_Sequence.new(gi[1]).fetch_sequence.to_s
  end
  puts fasta_contents
  puts "\n\n\n\n"
  File.open("../results/protein/partitions_pre_alignment/#{partition}_unaligned.fasta", 'w') { |f| f.write(fasta_contents) }
end

#align fasta files
Dir.foreach("../results/protein/partitions_pre_alignment/") do |file|
  if file.include? ".fasta"
    outfile = file.sub(/unaligned/, 'aligned')
    #platform-specific calls to muscle
    case
    when RUBY_PLATFORM.downcase.include?("linux")
      `./muscle3.8.31_i86linux64 -in ../results/protein/partitions_pre_alignment/#{file} -out ../results/protein/partitions_aligned/#{outfile}`
    when RUBY_PLATFORM.downcase.include?("darwin")
      `./muscle3.8.31_i86darwin64 -in ../results/protein/partitions_pre_alignment/#{file} -out ../results/protein/partitions_aligned/#{outfile}`
    end
  end
end
=end

################################################################################
# CREATE TNT INPUT AND BUILD TREES IN TNT
################################################################################
#keep total length of all partitions, defined outside of loop
total_nchar = 0

Dir.foreach("../results/protein/partitions_aligned") do |file|

  #ensure only fasta files are used in the directory - just as a precaution
  if file.include? ".fasta"

    #first, convert fasta files to tnt data matrices
    outfile = file.sub(/_aligned\.fasta/, '.tnt')
    input_file = "../results/protein/partitions_aligned/#{file}"
    output_file = "../results/protein/tnt_input/#{outfile}"

    #Instantiate strings to store in TNT files
    output = String.new()
    matrix = String.new()

    #Build TNT datafile
    output << "nstates prot;\n"
    output << "xread\n"

    input = Bio::FastaFormat.open(input_file)
    nchar = 0
    ntax = 0

    input.each do |entry|
      nchar = entry.seq.length() if nchar == 0
      ntax += 1
      taxname = Bio::GenPept.new(ncbi_fetch.sequence(entry.entry_id)).organism().tr!(" ", "_")
      matrix << "'#{taxname}' "
      matrix << entry.seq
      matrix << "\n"
    end
    output << "#{nchar} #{ntax}\n"
    output << matrix

    output << ";\n\n\nproc /;\ncomments 0\n;"

    puts output

    #File.open(output_file, "w") { |f| f.write(output) }

    #now, execute 'myscript.run' within tnt to build the trees
    case
    when RUBY_PLATFORM.downcase.include?("linux")
      `#{$tnt_path}/tnt mxram 2000 , p ../results/protein/tnt_input/#{outfile} , cd tnt , sect:slack 10 , myscript , zzz ,`
    when RUBY_PLATFORM.downcase.include?("darwin")
      `#{$tnt_path}/tnt.command mxram 2000 , p ../results/protein/tnt_input/#{outfile} , cd tnt , sect:slack 10 , myscript , zzz ,`
    end

    #save tree (parenthetical notation '.tre') and all tnt output ('.out') for current tree
    tnt_output = File.read("tnt/myscript.out")
    File.open("../results/protein/tnt_output/#{file.sub(/_aligned\.fasta/, '_tnt_output.out')}", 'w') {|f| f.write(tnt_output)}
    current_tree = File.read("tnt/myscript.tre")
    File.open("../results/protein/tnt_output/trees/#{file.sub(/_aligned\.fasta/, '_tree.tre')}", 'w') {|f| f.write(current_tree)}

    total_nchar += nchar
  end
end


################################################################################
# CREATE CONSENSUS TREE WITH PARTITIONED BREMER SUPPORT
################################################################################
#create interleaved TNT data matrix of all partitions
pbtnt = String.new()
int_matrix = String.new()

total_ntax = results.first[1].size

pbtnt << "nstates prot;\n"
pbtnt << "xread\n"
pbtnt << "#{total_nchar} #{total_ntax}\n\n"

results.each_key do |partition|
  current_matrix = String.new()
  current_matrix << "&[prot]\n"
  input = Bio::FastaFormat.open("../results/protein/partitions_aligned/#{partition}_aligned.fasta")
  input.each do |entry|
    taxname = Bio::GenPept.new(ncbi_fetch.sequence(entry.entry_id)).organism().tr!(" ", "_")
    current_matrix << "#{taxname} "
    current_matrix << entry.seq
    current_matrix << "\n"
  end
  int_matrix << current_matrix
  int_matrix << "\n"
end

pbtnt << int_matrix
pbtnt << ";\nproc /;"

File.open("../results/protein/tnt_input/pbsup.tnt", "w") { |f| f.write(pbtnt) }

#Calculate consensus tree with partitioned Bremer support
case
when RUBY_PLATFORM.downcase.include?("linux")
  `#{$tnt_path}/tnt cd tnt , mxram 2000 , p ../../results/protein/tnt_input/pbsup.tnt , pbsup #{results.size} , export * pbsup.nex , zzz ,`
when RUBY_PLATFORM.downcase.include?("darwin")
  `#{$tnt_path}/tnt.command cd tnt , mxram 2000 , p ../../results/protein/tnt_input/pbsup.tnt , pbsup #{results.size} , export * pbsup.nex , zzz ,`
end

#save consensus tree (parenthetical notation '.tre') and all tnt output ('.out')
tnt_output = File.read("tnt/pbs.out")
File.open("../results/protein/tnt_output/pbs.out", 'w') {|f| f.write(tnt_output)}
current_tree = File.read("tnt/pbs.tre")
File.open("../results/protein/tnt_output/trees/pbs.tre", 'w') {|f| f.write(current_tree)}

puts "Protein analysis completed - proceeding with Nucleotide analysis..."

################################################################################
# FINALLY, MAKE TREES FOR NUCLEOTIDE SEQUENCES...
################################################################################
=begin
################################################################################
# CREATE AND ALIGN FASTA FILES
################################################################################
# Create a fasta file for each partition, pulling sequences for all GIs from NCBI

nucleotide_results.each_pair do |partition, sequences|
  fasta_contents = String.new()
  sequences.each do |gi|
    pp gi
    if gi[1].scan(/^.{11}?join\(.+/).length > 0
      splice = Splice.new(gi[1], gi[0])#gi[0] is species, gi[1] is sequence...
      splice.parse_components()
      fasta_contents << splice.make_fasta().to_s
      puts "\nMade FASTA for transsplice gene #{gi}"
    else
      fasta_contents << Get_NCBI_Sequence.new(gi[1]).fetch_sequence.to_s
    end
  end
  puts fasta_contents
  puts "\n\n\n\n"
  File.open("../results/nucleotide/partitions_pre_alignment/#{partition}_unaligned.fasta", 'w') { |f| f.write(fasta_contents) }
end

#align fasta files
Dir.foreach("../results/nucleotide/partitions_pre_alignment/") do |file|
  if file.include? ".fasta"
    outfile = file.sub(/unaligned/, 'aligned')
    #platform-specific calls to muscle
    case
    when RUBY_PLATFORM.downcase.include?("linux")
      `./muscle3.8.31_i86linux64 -in ../results/nucleotide/partitions_pre_alignment/#{file} -out ../results/nucleotide/partitions_aligned/#{outfile}`
    when RUBY_PLATFORM.downcase.include?("darwin")
      `./muscle3.8.31_i86darwin64 -in ../results/nucleotide/partitions_pre_alignment/#{file} -out ../results/nucleotide/partitions_aligned/#{outfile}`
    end
  end
end
=end

=begin
################################################################################
# CREATE TNT INPUT AND BUILD TREES IN TNT
################################################################################
#keep total length of all partitions, defined outside of loop
total_nchar = 0

Dir.foreach("../results/nucleotide/partitions_aligned") do |file|

  #ensure only fasta files are used in the directory - just as a precaution
  if file.include? ".fasta"

    #first, convert fasta files to tnt data matrices
    outfile = file.sub(/_aligned\.fasta/, '.tnt')
    input_file = "../results/nucleotide/partitions_aligned/#{file}"
    output_file = "../results/nucleotide/tnt_input/#{outfile}"

    #########################################Skipping this block prevents rebuilding TNT data files
    #Instantiate strings to store in TNT files
    output = String.new()
    matrix = String.new()

    #Build TNT datafile
    output << "nstates dna;\n"
    output << "xread\n"

    input = Bio::FastaFormat.open(input_file)
    nchar = 0
    ntax = 0

    input.each do |entry|
      nchar = entry.seq.length() if nchar == 0
      ntax += 1
      if entry.entry_id.scan(/.*?gi.*?/).length > 0
        taxname = Bio::GenPept.new(ncbi_fetch.sequence(entry.entry_id)).organism().tr!(" ", "_")
      else
        taxname = entry.entry_id
      end
      taxname = Bio::GenPept.new(ncbi_fetch.sequence(entry.entry_id)).organism().tr!(" ", "_")
      matrix << "'#{taxname}' "
      matrix << entry.seq
      matrix << "\n"
    end
    output << "#{nchar} #{ntax}\n"
    output << matrix

    output << ";\n\n\nproc /;\ncomments 0\n;"

    puts output
    File.open(output_file, "w") { |f| f.write(output) }

    #now, execute 'myscript.run' within tnt to build the trees
    case
    when RUBY_PLATFORM.downcase.include?("linux")
      `#{$tnt_path}/tnt mxram 2000 , p ../results/nucleotide/tnt_input/#{outfile} , cd tnt , sect:slack 10 , myscript , zzz ,`
    when RUBY_PLATFORM.downcase.include?("darwin")
      `#{$tnt_path}/tnt.command mxram 2000 , p ../results/nucleotide/tnt_input/#{outfile} , cd tnt , sect:slack 10 , myscript , zzz ,`
    end

    #save tree (parenthetical notation '.tre') and all tnt output ('.out') for current tree
    tnt_output = File.read("tnt/myscript.out")
    File.open("../results/nucleotide/tnt_output/#{file.sub(/_aligned\.fasta/, '_tnt_output.out')}", 'w') {|f| f.write(tnt_output)}
    current_tree = File.read("tnt/myscript.tre")
    File.open("../results/nucleotide/tnt_output/trees/#{file.sub(/_aligned\.fasta/, '_tree.tre')}", 'w') {|f| f.write(current_tree)}
    total_nchar += nchar #UNCOMMENT FOR DEPLOYMENT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  end
end


################################################################################
# CREATE CONSENSUS TREE WITH PARTITIONED BREMER SUPPORT
################################################################################
#create interleaved TNT data matrix of all partitions
pbtnt = String.new()
int_matrix = String.new()

total_ntax = nucleotide_results.first[1].size

pbtnt << "nstates dna;\n"
pbtnt << "xread\n"
pbtnt << "#{total_nchar} #{total_ntax}\n\n"

nucleotide_results.each_key do |partition|
  current_matrix = String.new()
  current_matrix << "&[dna]\n"
  input = Bio::FastaFormat.open("../results/nucleotide/partitions_aligned/#{partition}_aligned.fasta")
  input.each do |entry|
    if entry.entry_id.scan(/.*?gi.*?/).length > 0
      taxname = Bio::GenPept.new(ncbi_fetch.sequence(entry.entry_id)).organism().tr!(" ", "_")
    else
      taxname = entry.entry_id
    end
    current_matrix << "#{taxname} "
    current_matrix << entry.seq
    current_matrix << "\n"
  end
  int_matrix << current_matrix
  int_matrix << "\n"
end

pbtnt << int_matrix
pbtnt << ";\nproc /;"

File.open("../results/nucleotide/tnt_input/pbsup.tnt", "w") { |f| f.write(pbtnt) }

#Calculate consensus tree with partitioned Bremer support
case
when RUBY_PLATFORM.downcase.include?("linux")
  `#{$tnt_path}/tnt cd tnt , mxram 2000 , p ../../results/nucleotide/tnt_input/pbsup.tnt , pbsup #{results.size} , export * pbsup.nex , zzz ,`
when RUBY_PLATFORM.downcase.include?("darwin")
  `#{$tnt_path}/tnt.command cd tnt , mxram 2000 , p ../../results/nucleotide/tnt_input/pbsup.tnt , pbsup #{results.size} , export * pbsup.nex , zzz ,`
end

#save consensus tree (parenthetical notation '.tre') and all tnt output ('.out')
tnt_output = File.read("tnt/pbs.out")
File.open("../results/nucleotide/tnt_output/pbs.out", 'w') {|f| f.write(tnt_output)}
current_tree = File.read("tnt/pbs.tre")
File.open("../results/nucleotide/tnt_output/trees/pbs.tre", 'w') {|f| f.write(current_tree)}

puts "Analysis completed"
=end
