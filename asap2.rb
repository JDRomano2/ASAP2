#!/usr/bin/env ruby

require "rubygems"
require "curses"
include Curses

$main_menu = {"Begin new analysis" => :new_session,
             "Continue current analysis" => :resume,
             "Load old analysis" => :load,
             "View software license" => :license,
             "Quit" => :quitProg}

def writeScreen()
  setpos(2,0)
  main_menu_string = ""
  $main_menu.keys.map.with_index { |x, i| main_menu_string << "   #{i+1}.  #{x}\n\n" }
  addstr(main_menu_string)
end

def quitProg
  close_screen
  exit
end

#### MAIN ####

init_screen
begin
  crmode
  noecho
  curs_set(0)
  writeScreen
  refresh
  sleep(2)
  quitProg
  ## MAIN contents go here
ensure
  close_screen
end

## END MAIN ##
