#!/usr/bin/env ruby

#NEXT: FIX INVALID STRINGS IN NUC_REFERENCE_SEQUENCES.JSON!

################################################################################
#
# asap2_6_cross_rf_values.rb
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

################################################################################
# DETERMINE ROBINSON-FOULDS DISTANCE BETWEEN CORRESPONDING PROT/NUC TREES
################################################################################

#NOTE!!!!!
#MAY NOT WORK, BECAUSE NUCLEOTIDE AND PROTEIN TREES ARE DERIVED FROM DIFFERENT
#INTERLEAVED DATA MATRICES!!!!

nuc_treefiles = Array.new()
Dir.foreach('../results/nucleotide/tnt_output/trees') {|treefile| nuc_treefiles.push(treefile) unless treefile.scan(/.*_tree\.tre.*/) == []}
prot_treefiles = Array.new()
Dir.foreach('../results/protein/tnt_output/trees') {|treefile| prot_treefiles.push(treefile) unless treefile.scan(/.*_tree\.tre.*/) == []}
rows = Array.new()
firstrow = ['Gene', 'RF Value']
rows.push(firstrow)
rfstring = String.new()

nuc_treefiles.each do |treefile1|
  current_row = ["#{treefile1}"]
  case
  when RUBY_PLATFORM.downcase.include?("linux")
    `#{$tnt_path}/tnt cd tnt , mxram 2000 , p ../../results/nucleotide/tnt_input/pbsup.tnt , p ../../results/nucleotide/tnt_output/trees/#{treefile1} , p ../../results/protein/tnt_output/trees/#{treefile1} , rfdistances 0 1 , zzz ,`
  when RUBY_PLATFORM.downcase.include?("darwin")
    `#{$tnt_path}/tnt.command cd tnt , mxram 2000 , p ../../results/nucleotide/tnt_input/pbsup.tnt , p ../../results/nucleotide/tnt_output/trees/#{treefile1} , p ../../results/protein/tnt_output/trees/#{treefile1} , rfdistances 0 1 , zzz ,`
  end
  rfvalue = IO.read('tnt/rflog.log').chomp
  transformed = (1 / (Math::E ** rfvalue.to_f))
  rfstring << "#{treefile1}   #{rfvalue}    #{transformed}\n"
  current_row.push(rfvalue)
  rows.push(current_row)
end

table = Terminal::Table.new :rows => rows
puts table
puts "\n\n"
puts "Gene name    RF value    Transformed RF value"
puts rfstring