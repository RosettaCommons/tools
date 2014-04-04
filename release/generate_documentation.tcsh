#!/bin/tcsh
#The purpose of this script is to generate the static wiki documentation for a weekly release.
#It runs from within prepare_weekly_release.bash
#I am ashamed that it is a tcsh script mixed in with bash scripts, but when Tim helpfully set up the environment for me he assumed I needed it set up in tcsh (my shell) instead of bash (used for the release for portability).
#author: Steven Lewis, smlewi@gmail.com
#intended to be run only on Contador, but can be safely edited for whatever other machine is used to create the weekly Rosetta release.

#this alias already exists in my .tcshrc
#alias rvm 'eval `~/bin/rvm.rb \!*`'

rvm use 1.9.3
gollum-site generate --base_path "./"