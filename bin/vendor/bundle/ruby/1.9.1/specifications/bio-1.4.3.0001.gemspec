# -*- encoding: utf-8 -*-

Gem::Specification.new do |s|
  s.name = "bio"
  s.version = "1.4.3.0001"

  s.required_rubygems_version = Gem::Requirement.new(">= 0") if s.respond_to? :required_rubygems_version=
  s.authors = ["BioRuby project"]
  s.date = "2013-05-24"
  s.description = "BioRuby is a library for bioinformatics (biology + information science)."
  s.email = "staff@bioruby.org"
  s.executables = ["bioruby", "br_biofetch.rb", "br_bioflat.rb", "br_biogetseq.rb", "br_pmfetch.rb"]
  s.extra_rdoc_files = ["KNOWN_ISSUES.rdoc", "README.rdoc", "README_DEV.rdoc", "RELEASE_NOTES.rdoc", "doc/Changes-1.3.rdoc", "doc/RELEASE_NOTES-1.4.0.rdoc", "doc/RELEASE_NOTES-1.4.1.rdoc", "doc/RELEASE_NOTES-1.4.2.rdoc"]
  s.files = ["bin/bioruby", "bin/br_biofetch.rb", "bin/br_bioflat.rb", "bin/br_biogetseq.rb", "bin/br_pmfetch.rb", "KNOWN_ISSUES.rdoc", "README.rdoc", "README_DEV.rdoc", "RELEASE_NOTES.rdoc", "doc/Changes-1.3.rdoc", "doc/RELEASE_NOTES-1.4.0.rdoc", "doc/RELEASE_NOTES-1.4.1.rdoc", "doc/RELEASE_NOTES-1.4.2.rdoc"]
  s.homepage = "http://bioruby.org/"
  s.rdoc_options = ["--main", "README.rdoc", "--title", "BioRuby API documentation", "--exclude", "\\.yaml\\z", "--line-numbers", "--inline-source"]
  s.require_paths = ["lib"]
  s.rubyforge_project = "bioruby"
  s.rubygems_version = "1.8.23"
  s.summary = "Bioinformatics library"

  if s.respond_to? :specification_version then
    s.specification_version = 3

    if Gem::Version.new(Gem::VERSION) >= Gem::Version.new('1.2.0') then
    else
    end
  else
  end
end
