#!/usr/bin/ruby

##############################################
#
# inParanoid_organisms.rb
#
# Copyright: Copyright (C) 2013 by Joe Romano
# Contact: jdromano2@gmail.com
# License: GNU General Public License
#
##############################################
#
# File Description:
# => finds overlap of all organisms in all inParanoid groups
#
##############################################

require 'rubygems'
require 'bio'
require 'pp'
require 'yaml'

all_files = Array.new()

Dir.glob('../../data/inParanoid/*.xml') do |gene|
  single_file = Array.new()
  handle = File.open(gene, 'r')
  handle.readlines.each do |line|
    match = line.scan(/speclong="(.+)"\/\>/)
    if match.length > 0
      single_file.push(match)
    end
  end
  all_files.push(single_file.uniq)
end

all_files.each { |e| pp e }

organisms = all_files.inject(:&)

puts organisms