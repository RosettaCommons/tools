#!/usr/bin/python

from sys import argv,stdout
import string
from math import sqrt
from parse_options import parse_options
from generate_constraints import generate_constraints

generate_constraints( argv, " CA ", " CA ", 12.0 )
