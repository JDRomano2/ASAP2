#!/usr/bin/env ruby

#NEXT: FIX INVALID STRINGS IN NUC_REFERENCE_SEQUENCES.JSON!

################################################################################
#
# asap2_2_find_nucs.rb
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
# FIND NUCLEOTIDE SEQUENCES ENCODING PROTEIN SEQUENCES
################################################################################

results = JSON.parse( IO.read('../results/logs/reference_sequences.json') )

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
