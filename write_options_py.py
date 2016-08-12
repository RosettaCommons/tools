#python
import argparse

print "This script takes a file containing options like so:\n\
type name description default\n\
Boolean PumpAction use \"the PumpAction mover\" true\n\
\n\
Everything after the third item is parsed as the description\n\
\n\
It returns option defs to paste into options_rosetta.py, like so:\n\
\n\
		Option('', '',\n\
			default = '',\n\
			desc = ""\n\
		),\n\
\n\
		Option('PumpAction', 'Boolean',\n\
			default = 'true',\n\
			desc = \"use the PumpAction mover\"\n\
		),\n\
###############paste below line into options_rosetta.py########"

parser = argparse.ArgumentParser(description="Write options_rosetta.py contents")

parser.add_argument('input', help='input options list', type=str)

args = parser.parse_args()

#args.input

#expected line types
#Boolean PumpAction "use the PumpAction mover" true
#0        1         2-x    last
for line in open(args.input, 'r'):
    if line[0] == "#":
        pass

    splitline = line.split()

    if len(splitline) < 3:
        print "Line is too short!"
        print line
        quit()

    #process into named vars
    type = splitline[0]
    name = splitline[1]
    default = splitline[-1]

    #process multiple word desc
    desc = ""
    for each in splitline[2:-1]:
        desc = desc + each + " "
    desc = desc.rstrip()

    #print results
    print "		Option('" + name + "', '" + type + "',"
    print "			default = '" + default + "',"
    print "			desc = '" + desc + "'"
