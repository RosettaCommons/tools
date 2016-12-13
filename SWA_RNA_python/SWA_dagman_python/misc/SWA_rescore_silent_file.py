#!/usr/bin/python

################################################################################
### import modules
################################################################################
import argparse
import subprocess as sp
import multiprocessing as mp
import os
from os.path import exists, basename, dirname, isdir, isfile, realpath
from glob import glob


################################################################################
### GLOBAL DEFINITIONS
################################################################################
DEFAULT_SILENT_FILE_NAME = 'region_*.out'
DEFAULT_OUTFILE = 'default.out'
DEFAULT_PATCH_SELECTORS = ['VIRTUAL_RIBOSE']
DEFAULT_WEIGHTS = 'stepwise/rna/rna_loop_hires_04092010.wts'
BACKUP_EXT = '.bak'
RNA_SCORE_OUT_EXT = '.rna_score'


################################################################################
### helper classes and functions
################################################################################
class Command(object):

    def __init__(self, executable):
        self.exe = executable
        self.args = []
        self.ifiles = []
        self.ofiles = []
        self.out = None
        self.err = None
        self.verbose = True

    def command_line(self, string = False):
        command = [self.exe] + self.args + self.ifiles + self.ofiles
        if string:
            command = ' '.join(command)
        return command

    def add_argument(self, argument, values = []):
        self.args.append(argument)
        if not isinstance(values, list):
            values = values.split(' ')
        self.args += [v for v in values]

    def add_infile(self, file):
        self.ifiles.append(file)

    def add_outfile(self, file):
        self.ofiles.append(file)

    def check_files(self):
        for file in self.ofiles:
            if not exists(file):
                if self.verbose:
                    print 'WARNING: %s does not exist !!!' % file
                return False
        for file in self.ifiles:
            if exists(file) and self.exe in ['mv', 'rm']:               
                if self.verbose:
                    print 'WARNING: %s still exists !!!' % file
                return False
        return True

    def check_error(self):
        if self.err and self.err != '':
            if self.verbose:
                print 'OUTPUT:', self.out
                print 'ERROR :', self.err
            return False
        return True

    def submit(self):
        if self.verbose:
            print 'Running command:', self.command_line(string=True)
        s = sp.Popen(self.command_line(), stdout=sp.PIPE, stderr=sp.PIPE)
        self.out, self.err = s.communicate()
        return ( self.check_error() and self.check_files() )


################################################################################
### rescore all silent files for a given target
################################################################################
def rescore(target, (filename, weights, patch_selectors, force, verbose)):
    
    ########################################################################
    ### change into target_dir and find all silent files
    ########################################################################
    working_dir = os.getcwd()
    os.chdir(target)    
    silent_files = sorted(glob(filename))
    
    ########################################################################
    ### iterate over all silent_files
    ########################################################################
    for idx, silent_file in enumerate(silent_files, start=1):

        silent_file_rna_score = silent_file + RNA_SCORE_OUT_EXT
        if verbose:
            print ( '[%d/%d] Rescoring silent file: %s' 
                    % (idx, len(silent_files), realpath(silent_file)) )

        ####################################################################
        ### if force, remove existing .rna_score file, otherwise, continue
        ####################################################################
        if exists(silent_file_rna_score):
            if not force:  continue
            command = Command('rm')
            command.add_infile(silent_file_rna_score)
            command.verbose = verbose
            if not command.submit():  continue

        ####################################################################
        ### remove DEFAULT_OUTFILE to avoid errors, only after backing up
        ####################################################################
        if exists( DEFAULT_OUTFILE ):
            command = Command('rm')
            command.add_infile(DEFAULT_OUTFILE)
            command.verbose = verbose
            if not command.submit():  continue

        ####################################################################
        ### run rna_score on silent_file
        ####################################################################
        command = Command('rna_score')
        command.add_argument('-score:weights', values = weights)
        command.add_argument('-patch_selectors', values = patch_selectors)
        command.add_argument('-silent', values = silent_file)
        command.verbose = verbose
        if not command.submit():  continue

        ####################################################################
        ### move default.out to *.out.rna_score
        ####################################################################
        command = Command('mv')
        command.add_infile(DEFAULT_OUTFILE)
        command.add_outfile(silent_file_rna_score)
        command.verbose = verbose
        if not command.submit():  continue

    ########################################################################
    ### change back to working directory
    ########################################################################
    os.chdir(working_dir)
    
    return True

################################################################################
### main script
################################################################################
if __name__=='__main__':


    ################################################################################
    ### parse arguments
    ################################################################################
    parser = argparse.ArgumentParser(
        description='Rescore all silent files for all targets in a given SWA run.'
    )

    parser.add_argument(
        'inpaths',
        nargs='+',
        help='List of paths to target runs.'
    )
    parser.add_argument(
        '-t','--targets',
        nargs='+',
        help='List of targets.',
        default=[]
    )
    parser.add_argument(
        '-s','--silent_file_name',
        help='Silent file to rescore.',
        default=DEFAULT_SILENT_FILE_NAME
    )
    parser.add_argument(
        '-weights',
        help='weights file to be used by rna_score.',
        default=DEFAULT_WEIGHTS
    )
    parser.add_argument(
        '-patch_selectors',
        nargs='*',
        help='additional patch_selectors that may be required.',
        default=DEFAULT_PATCH_SELECTORS
    )
    parser.add_argument(
        '-f','--force',
        help='Remove existing .rna_score files.',
        action='store_true'
    )
    parser.add_argument(
        '-j','--nproc',
        help='Number of jobs to run in parallel',
        default=(mp.cpu_count()-1)
    )

    args = parser.parse_args()
    
    inpaths = sorted( args.inpaths )
    user_targets = sorted( args.targets )
    silent_file_name = args.silent_file_name
    weights = args.weights
    patch_selectors = args.patch_selectors
    force = args.force
    nproc = int( args.nproc )

    ########################################################################
    ### checks and initializations 
    ########################################################################
    assert( all( exists( inpath ) for inpath in inpaths ) )
    working_dir = os.getcwd()

    ########################################################################
    ### change into inpaths and get targets
    ########################################################################
    for idx, inpath in enumerate(inpaths, start=1):
        
        os.chdir( inpath )
        print '\n[%d/%d] Rescoring targets in: %s' % (idx,len(inpaths),inpath)

        ########################################################################
        ### find all targets in inpath
        ########################################################################
        targets = sorted(filter(lambda x: isdir(x), glob('*')))
        if len(user_targets):
            targets = filter(lambda x: x in user_targets, targets)

        ########################################################################
        ### rescore all targets found in inpath
        ########################################################################
        if nproc > 1:
            
            ###################################################################
            ### Parallel version
            ###################################################################
            lock = mp.Lock()
            pool = mp.Pool(processes=nproc)
            
            arg = (silent_file_name, weights, patch_selectors, force, True)
            out = [pool.apply_async(rescore,args=(t,arg)) for t in targets]
            out = [o.get() for o in out]
            
            pool.close()
            pool.join()
           
        else:
            
            ###################################################################
            ### iterate over all target_dirs for a given inpath
            ###################################################################
            arg = (silent_file_name, weights, patch_selectors, force, True)
            out = [rescore(t,arg) for t in targets]           

        print out
        
        ########################################################################
        ### change back into working directory
        ########################################################################
        os.chdir( working_dir )

