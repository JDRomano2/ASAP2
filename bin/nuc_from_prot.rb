require 'bio'

puts "enter the GI: "
gi = gets


ncbi_fetch = Bio::NCBI::REST::EFetch.new()
Bio::NCBI.default_email = 'jdromano@uvm.edu'

record = Bio::GenPept.new(ncbi_fetch.sequence(gi))
CDS = record.features[(record.features.length - 1)]
codedby = nil
CDS.qualifiers.each { |a_qualifier| codedby = a_qualifier.value if a_qualifier.qualifier == "coded_by" }
if codedby == nil
  raise ArgumentError.new("Could not find nucleotide sequence that codes for #{gi}")
end
codedby.sub!(/\..*$/, '')
