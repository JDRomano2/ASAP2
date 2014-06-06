require 'rubygems'
require 'bio'

Bio::NCBI.default_email = "jdromano@uvm.edu"
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

test_gi = "NC_012920.1:3307..4262"
species = "Foo bar"

=begin
if test_gi.scan(/^(.+):<?(\d+)..(\d+)$/).length > 0
  splice = Splice.new(test_gi, species)#gi[0] is species, gi[1] is sequence...
  splice.parse_components()
  puts splice.make_fasta().to_s
  puts "\nMade FASTA for transsplice gene #{gi}"
=end
if test_gi.scan(/^(.+):<?(\d+)..(\d+)$/).length > 0
  seq = $1
  start = $2
  finish = $3
  
  splice_string = "#{start}..#{finish}"
  seq_all = Bio::Sequence.auto(ncbi_fetch.sequence(seq, "fasta").gsub!(/^>.+\n/, ''))
  fasta_definition = Get_NCBI_Sequence.new(seq).fetch_sequence.to_s.lines.first
  seq_trimmed = seq_all.splice(splice_string)

  #seq_all = Get_NCBI_Sequence.new(seq).fetch_sequence.to_s
  #fasta_definition = seq_all.lines.first
  #seq_trimmed = seq_all[(start.to_i)..(finish.to_i)]

  puts "it worked!"
  puts "\n\n"
  puts seq_all
  puts "\n\n"
  puts fasta_definition
  puts seq_trimmed
else
  puts "error!"
end