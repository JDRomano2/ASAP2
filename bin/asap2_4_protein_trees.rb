#!/usr/bin/env ruby

################################################################################
#
# asap2_4_protein_trees.rb
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
    rescue
      puts "There was an error getting the sequence.. retrying."
      retry
    end
    return @sequence
  end

  attr_accessor :sequence
end
results = JSON.parse( IO.read('../results/logs/reference_sequences.json') )

`mkdir -p ../results/protein/partitions_pre_alignment`
`mkdir -p ../results/protein/partitions_aligned`
`mkdir -p ../results/protein/tnt_input`
`mkdir -p ../results/protein/tnt_output/trees`


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
      `./muscle3.8.31_i86linux64 -in ../results/protein/partitions_pre_alignment/#{file} -out ../results/protein/partitions_aligned/unsorted_#{outfile}`
    when RUBY_PLATFORM.downcase.include?("darwin")
      `./muscle3.8.31_i86darwin64 -in ../results/protein/partitions_pre_alignment/#{file} -out ../results/protein/partitions_aligned/unsorted_#{outfile}`
    end
    `python stable.py ../results/protein/partitions_pre_alignment/#{file} ../results/protein/partitions_aligned/unsorted_#{outfile} > ../results/protein/partitions_aligned/#{outfile}`
    `rm ../results/protein/partitions_aligned/unsorted_#{outfile}`
  end
end

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
      begin
        taxname = Bio::GenPept.new(ncbi_fetch.sequence(entry.entry_id)).organism().tr!(" ", "_")
      rescue Exception
        puts "HTTP error... retrying."
        retry
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


##################
exit
##################

################################################################################
# CREATE CONSENSUS TREE WITH PARTITIONED BREMER SUPPORT
################################################################################
#create interleaved TNT data matrix of all partitions
pbtnt = String.new()
int_matrix = String.new()

total_nchar = 16198 ################ WARNING!!!!!!!! CHANGE THIS BEFORE RUNNING AGAIN!!!

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

File.open("../results/protein/tnt_input/pbsup.tnt", "w") { |f| f.write(pbtnt) }

#Calculate consensus tree with partitioned Bremer support
case
when RUBY_PLATFORM.downcase.include?("linux")
  `#{$tnt_path}/tnt cd tnt , mxram 2000 , p ../../results/protein/tnt_input/pbsup.tnt , pbsup #{results.size} , zzz ,`
when RUBY_PLATFORM.downcase.include?("darwin")
  `#{$tnt_path}/tnt.command cd tnt , mxram 2000 , p ../../results/protein/tnt_input/pbsup.tnt , pbsup #{results.size} , zzz ,`
end

#save consensus tree (parenthetical notation '.tre') and all tnt output ('.out')
tnt_output = File.read("tnt/pbs.out")
File.open("../results/protein/tnt_output/pbs.out", 'w') {|f| f.write(tnt_output)}
current_tree = File.read("tnt/pbs.tre")
File.open("../results/protein/tnt_output/trees/pbs.tre", 'w') {|f| f.write(current_tree)}

################################################################################
# DETERMINE ROBINSON-FOULDS DISTANCE BETWEEN EACH TREE (NOT PBS TREE)
################################################################################

treefiles = Array.new()
Dir.foreach('../results/protein/tnt_output/trees') {|treefile| treefiles.push(treefile) unless treefile.scan(/.*_tree\.tre.*/) == []}
rows = Array.new()
firstrow = ['']
treefiles.each {|treefile| firstrow.push(treefile)}
rows.push(firstrow)

rf_array = Array.new()
used = Array.new()

rfstring = String.new()

treefiles.each do |treefile1|
  current_row = ["#{treefile1}"]
  treefiles.each do |treefile2|
    if treefile1 == treefile2
      rfvalue = "n/a"
    else
      case
      when RUBY_PLATFORM.downcase.include?("linux")
        `#{$tnt_path}/tnt cd tnt , mxram 2000 , p ../../results/protein/tnt_input/pbsup.tnt , p ../../results/protein/tnt_output/trees/#{treefile1} , p ../../results/protein/tnt_output/trees/#{treefile2} , rfdistances 0 1 , zzz ,`
      when RUBY_PLATFORM.downcase.include?("darwin")
        `#{$tnt_path}/tnt.command cd tnt , mxram 2000 , p ../../results/protein/tnt_input/pbsup.tnt , p ../../results/protein/tnt_output/trees/#{treefile1} , p ../../results/protein/tnt_output/trees/#{treefile2} , rfdistances 0 1 , zzz ,`
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

File.open("../results/protein/rf_distances.txt", 'w') {|f| f.write(table)}
File.open("../results/protein/rf_distances.txt", 'a') {|f| f.write(rfstring)}

puts "Protein analysis completed"
