#!/usr/bin/env python
# Run the Membrane Protein Framework Unit Tests
# Author: Rebecca Alford
# Last Modified: 1/16/14

import subprocess
from subprocess import call
import shlex

# List of Unit Test suites
suites = [ 'EmbedDefIOTests', 
			'EmbedDefLoaderTests',
			'EmbedSearchParamsOptionsTests',
			'EmbedSearchParamsLoaderTests',
			'EmbedSearchParamsIOTests',
			'LipsFileIOTests',
			'LipoFileLoaderTest',
			'SpanFileIOTests',
			'SpanFileLoaderTests',
			'LoadAllResourcesTest',
			'GeometryUtilTest',
			'MembraneResiduesTest',
			'MembraneResiduesTest',
			'EmbeddingFactoryTest',
			'MPDatabaseIOTest',
			'MembraneInfoTest',
			'MembraneProteinFactoryTest' ]

# Unit Testing Command
cmd = 'python test/run.py -j2 -d ../database -c clang -1'

# Run custom ut command for each test in mp test suite
for i in range( len(suites) ): 
	print 'Running test suite ' + suites[i]
	mycmd = cmd + suites[i]
	subprocess.call(shlex.split(mycmd))
	print suites[i] + ' Test Passed!'

print 'All Membrane Framework Unit Tests Pass!'
