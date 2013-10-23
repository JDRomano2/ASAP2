#!/usr/bin/env ruby

#Class accepts a CDS from a transspliced gene (i.e., protein assembled by many nuc sequences across genome) and splices it into correct nuc. sequence
#
#Modeled off of following commands:
#
#seq = Bio::Sequence.auto(ncbi_fetch.sequence("KB031078.1", format = 'fasta'))
#puts seq.splice('complement(join(1185878..1186086))')
#
#in other words, create a Bio::Sequence object, store the FULL CHROMOSOMAL ASSEMBLY IN IT, then splice ranges of sequences in next command...

require 'rubygems'
require 'bio'
require 'json'
require 'terminal-table'
require 'pp'

class Splice
  attr_accessor :fasta
  attr_accessor :accession
  attr_accessor :seq_spliced

  def initialize(complement, source)
    @ncbi_fetch = Bio::NCBI::REST::EFetch.new()
    @complement = complement
    @source = source
    @organism = Bio::GenPept.new(ncbi_fetch.sequence(@source)).organism().tr!(" ", "_")
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
