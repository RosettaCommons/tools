#!/usr/bin/env python

from sys import argv, stdout

def Help():
    print argv[0], ' <type> <variable name 1 > <variable name 2> ... '
    exit( 0 )

if len( argv ) <= 1: Help()

var_type = argv[1]
if var_type == 'Size': var_type = 'core::Size'
if var_type == 'Real': var_type = 'core::Real'
if var_type == 'string': var_type = 'std::string'

default_vals = {'bool':'false','core::Size':'0','core::Real':'0.0', 'std::string':'""'}
#if var_type not in default_vals.keys():    Help()

vars = argv[2:]

print "Put following after public: in .hh"
print
for var in vars:
    print "\t\tvoid set_%s( %s const & setting ){ %s_ = setting; }" % (var, var_type, var )
    print "\t\t%s %s() const { return %s_; }" % (var_type, var, var )
    print
print
print
print "Put following after private: in .hh!"
print
#stdout.write( "\t%s %s" % (var_type, vars[0]),)
#for var in vars[1:]: stdout.write( ", %s_" % var )
#stdout.write( ";\n" )
for var in vars:
    print "\t\t%s %s_;" % (var_type, var )
print
print
print
print "Put following after class_definition"
print
for var in vars:
    if var_type in default_vals.keys():
        print '\t\t%s_( %s ),' % (var, default_vals[ var_type ] )


