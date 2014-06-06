#!/usr/bin/env ruby

################################################################################
#
# asap2_5_nuc_trees.rb
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
    begin
      @sequence = Bio::FastaFormat.new(ncbi_fetch.sequence(@accessor, "fasta"))
    rescue Exception
      puts "HTTP error... retrying."
      retry
    end
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
    begin
      @seq = Bio::Sequence.auto(@ncbi_fetch.sequence(@accession, "fasta").gsub!(/^>.+\n/, '')) #gsub removes header line from fasta
    rescue Exception
      puts "HTTP error... retrying."
      retry
    end
    @range_string = @ranges.join(",")
    if (@complement.scan(/^.*.{11}join\(([0-9A-Z]+)\..*$/).length > 0)
      @seq_spliced = @seq.splice("complement(join(#{@range_string}))")
    else
      @seq_spliced = @seq.splice("join(#{@range_string})")
    end
    @fasta = Bio::FastaFormat.new(">#{@organism} transsplice mRNA coding protein #{@accession}\n#{@seq_spliced}")
  end
end

nucleotide_results = JSON.parse( IO.read('../results/logs/nuc_reference_sequences.json'))

`mkdir -p ../results/nucleotide/partitions_pre_alignment`
`mkdir -p ../results/nucleotide/partitions_aligned`
`mkdir -p ../results/nucleotide/tnt_input`
`mkdir -p ../results/nucleotide/tnt_output/trees`


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
    elsif gi[1].scan(/^complement\((\w+.\d+):<?>?(\d+)..<?>?(\d+)\)/).length > 0
      seq = $1
      start = $2
      finish = $3

      splice_string = "#{start}..#{finish}"
      begin
        seq_all = Bio::Sequence.auto(ncbi_fetch.sequence(seq, "fasta").gsub!(/^>.+\n/, ''))
      rescue Exception
        puts "HTTP error... retrying."
        retry
      end
      fasta_definition = Get_NCBI_Sequence.new(seq).fetch_sequence.to_s.lines.first
      seq_trimmed = seq_all.complement.splice(splice_string)
      fasta_contents << fasta_definition
      fasta_contents << seq_trimmed
      fasta_contents << "\n"
    elsif gi[1].scan(/^(.+):<?>?(\d+)..<?>?(\d+)$/).length > 0
      seq = $1
      start = $2
      finish = $3
  
      splice_string = "#{start}..#{finish}"
      begin
        seq_all = Bio::Sequence.auto(ncbi_fetch.sequence(seq, "fasta").gsub!(/^>.+\n/, ''))
      rescue Exception
        puts "HTTP error... retrying."
        retry
      end
      fasta_definition = Get_NCBI_Sequence.new(seq).fetch_sequence.to_s.lines.first
      seq_trimmed = seq_all.splice(splice_string)
      fasta_contents << fasta_definition
      fasta_contents << seq_trimmed
      fasta_contents << "\n"
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
      `./muscle3.8.31_i86linux64 -in ../results/nucleotide/partitions_pre_alignment/#{file} -out ../results/nucleotide/partitions_aligned/unsorted_#{outfile}`
    when RUBY_PLATFORM.downcase.include?("darwin")
      `./muscle3.8.31_i86darwin64 -in ../results/nucleotide/partitions_pre_alignment/#{file} -out ../results/nucleotide/partitions_aligned/unsorted_#{outfile}`
    end
    `python stable.py ../results/nucleotide/partitions_pre_alignment/#{file} ../results/nucleotide/partitions_aligned/unsorted_#{outfile} > ../results/nucleotide/partitions_aligned/#{outfile}`
    `rm ../results/nucleotide/partitions_aligned/unsorted_#{outfile}`
  end
end


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
        begin
          taxname = Bio::GenPept.new(ncbi_fetch.sequence(entry.entry_id)).organism().tr!(" ", "_")
        rescue Exception
          puts "HTTP error... retrying."
          retry
        end
      else
        taxname = entry.entry_id
      end
      if taxname.length < 3          ## Value is arbitrary... just very small
        puts "Error - No Species name provided for #{file}"
        puts "Enter species name manually:"
        taxname = gets.chomp()
      end
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

#############
exit 
#############

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

    # ONLY NEEDED FOR MTDNA + AD GENES!
    if taxname == "Otolemur_crassicaudatus"
      taxname = "Otolemur_garnettii"
    end
    if taxname == "Odobenus_rosmarus_rosmarus"
      taxname = "Odobenus_rosmarus_divergens"
    end
    if taxname == "Tupaia_belangeri"
      taxname = "Tupaia_chinensis"
    end
    if taxname == "Ceratotherium_simum"
      taxname = "Ceratotherium_simum_simum"
    end
    if taxname == "Trichechus_manatus"
      taxname = "Trichechus_manatus_latirostris"
    end
    if taxname == "Mustela_putorius"
      taxname = "Mustela_putorius_furo"
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
  `#{$tnt_path}/tnt cd tnt , mxram 2000 , p ../../results/nucleotide/tnt_input/pbsup.tnt , pbsup #{nucleotide_results.size} , export * pbsup.nex , zzz ,`
when RUBY_PLATFORM.downcase.include?("darwin")
  `#{$tnt_path}/tnt.command cd tnt , mxram 2000 , p ../../results/nucleotide/tnt_input/pbsup.tnt , pbsup #{nucleotide_results.size} , export * pbsup.nex , zzz ,`
end

#save consensus tree (parenthetical notation '.tre') and all tnt output ('.out')
tnt_output = File.read("tnt/pbs.out")
File.open("../results/nucleotide/tnt_output/pbs.out", 'w') {|f| f.write(tnt_output)}
current_tree = File.read("tnt/pbs.tre")
File.open("../results/nucleotide/tnt_output/trees/pbs.tre", 'w') {|f| f.write(current_tree)}

################################################################################
# DETERMINE ROBINSON-FOULDS DISTANCE BETWEEN EACH TREE (NOT PBS TREE)
################################################################################

treefiles = Array.new()
Dir.foreach('../results/nucleotide/tnt_output/trees') {|treefile| treefiles.push(treefile) unless treefile.scan(/.*_tree\.tre.*/) == []}
rows = Array.new()
firstrow = ['']
treefiles.each {|treefile| firstrow.push(treefile)}
rows.push(firstrow)
rfstring = String.new()

treefiles.each do |treefile1|
  current_row = ["#{treefile1}"]
  treefiles.each do |treefile2|
    if treefile1 == treefile2
      rfvalue = "n/a"
    else
      case
      when RUBY_PLATFORM.downcase.include?("linux")
        `#{$tnt_path}/tnt cd tnt , mxram 2000 , p ../../results/nucleotide/tnt_input/pbsup.tnt , p ../../results/nucleotide/tnt_output/trees/#{treefile1} , p ../../results/nucleotide/tnt_output/trees/#{treefile2} , rfdistances 0 1 , zzz ,`
      when RUBY_PLATFORM.downcase.include?("darwin")
        `#{$tnt_path}/tnt.command cd tnt , mxram 2000 , p ../../results/nucleotide/tnt_input/pbsup.tnt , p ../../results/nucleotide/tnt_output/trees/#{treefile1} , p ../../results/nucleotide/tnt_output/trees/#{treefile2} , rfdistances 0 1 , zzz ,`
      end
      rfvalue = IO.read('tnt/rflog.log').chomp
      transformed = (1 / (Math::E ** rfvalue.to_f))
      rfstring << "#{treefile1}   #{treefile2}   #{rfvalue}    #{transformed}\n"
    end
    current_row.push(rfvalue)
  end
  rows.push(current_row)
end

table = Terminal::Table.new :rows => rows
puts table
puts "\n\n"
puts "Gene 1    Gene 2    RF value    Transformed RF value"
puts rfstring

File.open("../results/nucleotide/rf_distances_nuc.txt", 'w') {|f| f.write(table)}
File.open("../results/nucleotide/rf_distances_nuc.txt", 'a') {|f| f.puts(rfstring)}

puts "Nucleotide analysis completed"
