#!/usr/bin/env python

from os import listdir
from os.path import exists, isfile, dirname, basename, expandvars
import string


def swa_rna_build_dag_parse():

    JOBS = {}
    SCRIPTS_PRE = {}
    SCRIPTS_POST = {}

    RNA_BUILD_DAG = 'rna_build.dag'

    if exists( RNA_BUILD_DAG ):
        rna_build_dag = open( RNA_BUILD_DAG, 'r' )
        flines = rna_build_dag.readlines()

        for line in flines:
            line = line.strip('\n').split(' ')
            if line[0] == 'JOB':
                JOBS[ line[1] ] = line[2]
            elif line[0] == 'SCRIPT':
                if line[1] == 'PRE':
                    SCRIPTS_PRE[ line[2] ] = ' '.join(line[3:])
                elif line[1] == 'POST':
                    SCRIPTS_POST[ line[2] ] = ' '.join(line[3:])
            elif line[0] == 'PARENT':
                continue
            else:
                continue

        rna_build_dag.close()


       # for k,v in JOBS.iteritems():
       #     print 'JOB: ', k
       #     print '>>>  ', v

       # for k,v in SCRIPTS_PRE.iteritems():
       #     print 'SCRIPT PRE: ', k
       #     print '>>>  ', v

       # for k,v in SCRIPTS_POST.iteritems():
       #     print 'SCRIPT POST: ', k
       #     print '>>>  ', v

    else:
        print "ERROR: ", RNA_BUILD_DAG, " does not exist!!!"

    return JOBS, SCRIPTS_PRE, SCRIPTS_POST
