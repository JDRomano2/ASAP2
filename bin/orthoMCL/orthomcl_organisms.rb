#!/usr/bin/ruby

##############################################
#
# orthomcl_organisms.rb
#
# Copyright: Copyright (C) 2013 by Joe Romano
# Contact: jdromano2@gmail.com
# License: GNU General Public License
#
##############################################
#
# File Description:
# => Finds overlap between all OrthoMCL .txt files in /data/OrthoMCL/
#
##############################################

require 'rubygems'
require 'bio'
require 'pp'
require 'yaml'

all_files = Array.new()

Dir.glob('../../data/OrthoMCL/*.txt') do |gene|
  single_file = Array.new()
  handle = File.open(gene, 'r')
  handle.readlines.each do |line|
    match = line.scan(/^\>.+\[(.+)\]/)
    if match.length > 0
      single_file.push(match)
    end
  end
  all_files.push(single_file.uniq)
end

all_files.each { |e| pp e }

organisms = all_files.inject(:&)

puts organisms
