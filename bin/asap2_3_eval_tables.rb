#!/usr/bin/env ruby

#NEXT: FIX INVALID STRINGS IN NUC_REFERENCE_SEQUENCES.JSON!

################################################################################
#
# asap2_3_eval_tables.rb
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

################################################################################
# CREATE TERMINAL TABLES FOR IDENTIFIED PROTEIN SEQUENCES
################################################################################

results = JSON.parse( IO.read('../results/logs/reference_sequences.json') )

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
      current_blast = `blastp -query blast_query.fasta -subject blast_subject.fasta -outfmt '6 evalue'`
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