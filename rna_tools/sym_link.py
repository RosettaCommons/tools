import os, errno
from os import chdir,getcwd
from os.path import basename,dirname
import glob
from sys import argv

#Might as well do a cross-check here that Rosetta exists.

#try:
#    rosetta_folder = os.environ["ROSETTA"]
#except KeyError:
#    print "You need to define ROSETTA as an environment variable, e.g., put "
#    print " export ROSETTA='/Users/rhiju/src/rosetta/' "
#    print "in your .bashrc or .bash_profile script."

#if not os.access('/path/to/folder', os.W_OK): exit( 0 )

CWD = getcwd()
chdir( dirname( argv[0] ) )

for f in glob.glob('./bin/*'):
    os.remove(f)

for dirpath, dirnames, filenames in os.walk('./'):
    if len( basename(dirpath) ) == 0 or basename(dirpath)[0] == '.' or basename(dirpath) == 'bin':
        continue
    for f in filenames:
        filename = os.path.join(dirpath, f)
        if filename[-3:] == '.py':
            #print filename
            file1 = '../' + filename
            file2 = './bin/' + f
            try:
                os.symlink(file1, file2)
            except OSError as e:
                if e.errno == errno.EEXIST:
                    os.remove(file2)
                    os.symlink(file1, file2)

chdir( CWD )
