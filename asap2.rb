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

def getUserChoice()
  selection = nil
  choice = getch
  case choice
  when '1'
    selection = :new_session
    addstr("You selected '1'!")
    sleep(0.5)
  when '2'
    selection = :resume
  when '3'
    selection = :load
  when '4'
    selection = :license
  when '5'
    selection = :quitProg
  end
  selection
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
  next_action = getUserChoice
  quitProg
  ## MAIN contents go here
ensure
  close_screen
end

## END MAIN ##
