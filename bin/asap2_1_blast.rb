#!/usr/bin/env ruby

#NEXT: FIX INVALID STRINGS IN NUC_REFERENCE_SEQUENCES.JSON!

################################################################################
#
# asap2_1_blast.rb
#
# Copyright: Copyright (C) 2013 by Joe Romano
# Contact: jdromano2@gmail.com
# License: GNU General Public License
#
################################################################################
#
# File Description:
# => ASAP2 - Automated Simultaneous Analysis Phylogenetics; Ver. 2
# => Accept query GIs and run BLAST analyses to find putative orthologues
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