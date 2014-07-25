#!/usr/bin/env ruby

################################################################################
#
# asap2.rb
#
# Copyright: Copyright (C) 2014 by Joe Romano
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
require 'bio'
require 'optparse'
require 'optparse/time'
require 'ostruct'
require 'fileutils' # For creating missing directories
require 'json'

###########################################
# INITIALIZE VARIABLES                    #
###########################################

$TNT_DIR = "/usr/local/bin"

###########################################
# DEFINE FUNCTIONS                        #
###########################################

def colorize(text, color_code)
  # for ANSI escape codes, see `http://en.wikipedia.org/wiki/ANSI_escape_code'
  "\e[#{color_code}m#{text}\e[0m"
end

def red(text); colorize(text, 31); end
def green(text); colorize(text, 32); end
def blue(text); colorize(text, 34); end
def magenta(text); colorize(text, 35); end
# e.g., `puts blue('this is in blue')`

$str = nil # store user responses to prompts - global scope

def waitForKey
  $str = nil # just to be safe
  begin
    system("stty raw -echo")
    $str = STDIN.getc
  ensure
    system("stty -raw echo")
  end
end

###########################################
# PARSE OPTIONS                           #
###########################################

class OptparseJDR
  def self.parse(args)
    options = OpenStruct.new

    option_parser = OptionParser.new do |opts|
      opts.banner = "Usage: wrapper.rb [options]"

      opts.on("-h", "--help", "Display help message") do |h|
        puts opts
        exit
      end

      opts.on("-v", "--version", "Show version") do |v|
        puts "ASAP2: Version 1.0"
        exit
      end
    end

    option_parser.parse!(args)
    options
  end
end
options = OptparseJDR.parse(ARGV)

# Detect OS and proceed accordingly
require 'rbconfig'
def os
  @os ||= (
    host_os = RbConfig::CONFIG['host_os']
    case host_os
    when /mswin|msys|mingw|cygwin|bccwin|wince|emc/
      ':windows'
    when /darwin|mac os/
      ':macosx'
    when /linux/
      ':linux'
    when /solaris|bsd/
      ':unix'
    else
      raise RuntimeError, "unknown os: #{host_os.inspect}"
    end
  )
end

###########################################
# DRAW APPLICATION TO SCREEN              #
###########################################

welcome_message = <<-END
ASAP2 - Automated Simultaneous Analysis Phylogenetics, Ver. 2
written by Joe Romano
License: GNU General Public License

Full documentation may be found in README.md or online at:
  http://github.com/JDRomano2/ASAP2

If this is your first time running this application, you will be prompted to
enter local configuration options. You may edit the configuration at any time
by manually altering the configuration file at <DIRECTORY HERE>, or by
selecting the appropriate option at the main menu.

Please press any key to continue...
END

$main_menu_string = <<-END
  --MAIN MENU--

(1)   Define new set of query sequences (for a new analysis)

(2)   Find potential orthologues to query sequences

(3)   Find reference sequences for organisms that are in all partitions

(4)   Generate nucleotide reference sequences coding for protein reference sequences

(5)   Run MP phylogenetic analysis on reference sequences

(6)   Find Robinson-Foulds distances between trees

(c)   Edit configuration options

(q)   Quit ASAP2


END

# Clear the screen


def main_menu()
  system "clear" or system "cls"
  puts $main_menu_string

  #TODO: program should always store the name of the active analysis

  puts green("Current analysis: ") + red("{ #{$ACTIVE_ANALYSIS} }")

  waitForKey

  # respond to previously entered key
  case $str
  when '1'
    initializeQueries
  when '2'
    findOrthologues
  when '3'
    findRefSeqs
  when '4'
    nucsFromProts
  when '5'
    findTrees
  when '6'
    findRFs
  when 'c'
    doConfig
  when 'q'
    exit
  else
  end

end

def initializeQueries
  puts red("test")
  sleep 2
end

def doConfig
  FileUtils::mkdir_p '../config'
  system "clear" or system "cls"
  config_hash = {}

  puts "Your email address?"
  config_hash[:email] = gets.chomp
  
  puts "Location of TNT executable? (default is `/usr/bin/tnt`)"
  config_hash[:tnt] = gets.chomp
  config_hash[:tnt] = "/usr/bin/tnt" if (config_hash[:tnt] == "")

  File.open('../config/config.txt', 'w') { |f| f.write(config_hash.to_json) }
end

def loadConfig
  config_hash = JSON.parse(IO.read('../config/config.txt'))
  Bio::NCBI.default_email = config_hash[:email]
end

def loadAnalysis
  analysisValid == false
  until analysisValid
    puts "Please enter the name of the folder containing the desired analysis, or type 'new' for a new analysis"
    user_input = gets.chomp
    if File.file?('../data/#{user_input}')
      analysisValid = true
      # load the file
    elsif user_input == "new"
      analysisValid = true
      # make new analysis file and initialize values
    else
      puts "Analysis name not recognized - please try again"
    end
  end
end

###########################################
# FILE OPERATIONS / ADTs                  #
###########################################

# This part of the program contains functions that interact with data files,
# as well as the abstract data types that store/operate on the data.
# Each subsection of ASAP2 requires input and output - they are dealt with here.

class Tree
  def initialize(name)
    @name = name
    @tree = ""
  end

  attr_accessor :name
  attr_accessor :tree
end

class RF
  def initialize(tree1, tree2)
    @tree1 = tree1
    @tree2 = tree2
  end
end

class Analysis
  def initialize(name, is_new)
    @name = name
    @query_gis = []
    @orthologues = []
    @prot_ref_seqs = []
    @nuc_ref_seqs = []
    @trees = []
    @rfs = []
    if is_new
      #build file structure for new analysis
      #TODO
    end
  end

  attr_accessor :name
  attr_accessor :query_gis
  attr_accessor :orthologues
  attr_accessor :prot_ref_seqs
  attr_accessor :nuc_ref_seqs
  attr_accessor :trees

  def exportToJSON
  end

  def saveToBinary
  end

  def parseQueries
  end

  def findOrthologues
  end

  def findProtRefSeqs
  end

  def findNucRefSeqs
  end

  def buildTrees
  end

  def findRFs
  end

end

#TODO
def loadFromJSON
end

def loadFromBinary
end

###########################################
# RUN PROGRAM                             #
###########################################

if __FILE__ == $0

  $ACTIVE_ANALYSIS = nil

  if File.file?('../config/config.txt')
    loadConfig
  else
    puts red("\n\n\tRunning for the first time,\n\tproceeding to configuration")
    sleep 2
    doConfig
  end
  
  system "clear" or system "cls"

  puts welcome_message

  waitForKey

  os_symbol = os()
  if os_symbol == :windows
    at_exit { puts red("Windows is not currently supported - please try running again on a UNIX-based computer.") }
    exit
  end

  at_exit { system "clear" or system "cls" }

  main_menu() while true

end
