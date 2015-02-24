#!/usr/bin/env python
# :noTabs=true:

# This file is made available under the Rosetta Commons license.
# See http://www.rosettacommons.org/license
# (C) 199x-2007 University of Washington
# (C) 199x-2007 University of California Santa Cruz
# (C) 199x-2007 University of California San Francisco
# (C) 199x-2007 Johns Hopkins University
# (C) 199x-2007 University of North Carolina, Chapel Hill
# (C) 199x-2007 Vanderbilt University

## @file   DeployPyRosetta.py
## @brief  Setup pyrosetta build environment.
## @author Sergey Lyskov

import os, sys, time, platform, commands, stat, subprocess, shutil

from optparse import OptionParser, IndentedHelpFormatter

# Create global 'Platform' that will hold info of current system
if sys.platform.startswith("linux"): Platform = "linux" # can be linux1, linux2, etc
elif sys.platform == "darwin" : Platform = "macos"
elif sys.platform == "cygwin" : Platform = "cygwin"
else: Platform = "_unknown_"
PlatformBits = platform.architecture()[0][:2]



def execute(message, commandline, return_=False):
    print message
    print commandline
    (res, output) = commands.getstatusoutput('bash -c "%s"' % commandline)
    print output
    if res:
        print "\nEncounter error while executing: ", commandline
        if return_: return True
        else: sys.exit(1)

    return False


def main(args):
    ''' Script to Setup build environment for PyRosetta.
    '''
    parser = OptionParser(usage="usage: %prog [OPTIONS]")
    parser.set_description(main.__doc__)

    parser.add_option('--prefix',
      #default=os.path.join( os.path.expanduser("~"), "PyRosetta.Develop." + PlatformBits),
      default=os.path.abspath('./PyRosetta.Develop.' + PlatformBits),
      action="store",
      help="Path to where PyRosetta build environmnet should be deployed.",
    )

    parser.add_option("--rosetta-source",
      default='',  # os.path.abspath('./rosetta_source'),
      help="Path to Rosetta source code where bash script to build PyRosetta will be created.",
    )

    parser.add_option("-j", "--jobs",
      default=1,
      type="int",
      help="Number of processors to use on when building(default: 1)",
    )

    parser.add_option("--omit-cmake",
      default=False,
      action="store_true",
      help="Disable installation of cmake. Use if it is already installed.",
    )

    parser.add_option("--compiler",
      default='gcc',
      action="store",
      help="Default compiler that will be used to build PyRosetta. Default is 'gcc'.",
    )


    parser.add_option("--loader",
      default='curl -OL',
      action="store",
      help="Default name of wget/curl like loader programm with extra options if nessesary. Default is 'curl -LO'.",
    )


    parser.add_option("-d",
      action="store_true", dest="debug",
      help="Enable DEBUG mode. Disable 30 seconds wait before start working.",
    )

    parser.add_option("--debug",
      default=False,
      )


    (options, args) = parser.parse_args(args=args[1:])
    global Options;  Options = options

    # determing gcc version...
    (res, output) = commands.getstatusoutput('gcc -v')
    res = output.rfind('version 4.9')
    if res >= 0:
        print GCC_DowngradeMessage % output
        sys.exit(1)



    print 'Starting set-up of PyRosetta build environment in 30 seconds. Please inspect the settings and press ctrl-c if you want to abort...'

    prefix = os.path.abspath(Options.prefix)
    #data_dir = os.path.abspath('./data')
    working_dir = os.path.abspath('./build')
    source_dir = working_dir + '/source'
    rosetta_source_dir = os.path.abspath(Options.rosetta_source)

    i_cmake = 'cmake-2.8.12.2'  #'cmake-2.8.7'
    #i_pygccxml = 'pygccxml-1.0.0'
    #i_pyplusplus = 'Py++-1.0.0'
    i_boost_version = (1, 55, 0)  # (1, 54, 0)
    i_boost = 'boost_%s_%s_%s' % i_boost_version
    i_python_lib = 'python%s.%s' % sys.version_info[:2]
    i_BuildPyrosetta = prefix +'/BuildPyRosetta.sh'

    print 'Destination path:', prefix
    #print 'Data dir:', data_dir
    print 'Working dir:', working_dir
    print 'Source dir: ', source_dir
    print 'Python dynamic lib to be used:', i_python_lib
    print 'PyRosetta building dir:', Options.rosetta_source if Options.rosetta_source else '[SKIPPED]'
    print 'Omit CMake install:', Options.omit_cmake
    print 'Compiler to use when buiilding PyRosetta:', Options.compiler
    print 'Debug:', Options.debug

    if not Options.debug : time.sleep(30)

    # Creating Dirs and set environment scripts.
    cd_work_env = 'source %s/PyRosettaBuildEnvironment.sh && cd %s && ' % (prefix, working_dir)

    if not os.path.isdir(prefix): os.makedirs(prefix)
    execute('Cleaning up working dir...', 'rm -rf %s' % working_dir);  os.makedirs(working_dir);  os.makedirs(source_dir)

    source = '#!/bin/sh\nexport PATH=%s/bin:$PATH\nexport PYTHONPATH=%s/lib/%s/site-packages:$PYTHONPATH\nexport LD_LIBRARY_PATH=%s/lib:$LD_LIBRARY_PATH\n' % (prefix, prefix, i_python_lib, prefix)
    file(prefix+'/PyRosettaBuildEnvironment.sh', 'w').write(source)
    execute('Creating environment script...', 'chmod +x %s/PyRosettaBuildEnvironment.sh' % prefix)

    execute('Testing source file...', cd_work_env + 'ls')

    if not Options.omit_cmake:
        execute('Downloading CMake...', 'cd %s && %s "http://www.cmake.org/files/v2.8/%s.tar.gz"' % (source_dir, Options.loader, i_cmake) )
        execute('Unpacking CMake...', 'cd %s && tar -vzxf %s/%s.tar.gz' % (working_dir, source_dir, i_cmake) )
        execute('Installing CMake...', 'cd %s/%s && ./configure --prefix=%s && make -j%s && make install' % (working_dir, i_cmake, prefix, Options.jobs) )

    # CVS version, now depricated
    # execute('Getting GCCXML...', 'cd %s && mkdir gccxml-cvs && cd gccxml-cvs && cvs -d :pserver:anoncvs@www.gccxml.org:/cvsroot/GCC_XML co -D 2012-10-20 gccxml' % (working_dir, ) )
    # execute('Creating GCCXML build dir...', cd_work_env + 'mkdir gccxml-cvs/gccxml-build')
    # execute('Configuring GCCXML...', cd_work_env + 'cd gccxml-cvs/gccxml-build && cmake ../gccxml -DCMAKE_INSTALL_PREFIX:PATH=%s' % prefix)
    # execute('Building GCCXML...', cd_work_env + 'cd gccxml-cvs/gccxml-build && make -j%s && make install' % Options.jobs )

    execute('Getting GCCXML...', 'cd %s && mkdir gccxml-git && cd gccxml-git && git clone https://github.com/gccxml/gccxml.git' % (working_dir, ) )  # git reset --hard 1bfbbe93ae16e97b010fb111305d169e3dcfd5a4
    execute('Creating GCCXML build dir...', cd_work_env + 'mkdir gccxml-git/gccxml-build')
    execute('Configuring GCCXML...', cd_work_env + 'cd gccxml-git/gccxml-build && cmake ../gccxml -DCMAKE_INSTALL_PREFIX:PATH=%s' % prefix)
    execute('Building GCCXML...', cd_work_env + 'cd gccxml-git/gccxml-build && make -j%s && make install' % Options.jobs )

    #execute('Unpacking PyGCCXML...', 'cd %s && unzip  %s/%s.zip' % (working_dir, data_dir, i_pygccxml) )
    #execute('Installing PyGCCXML...', cd_work_env + 'cd %s && python setup.py install --prefix=%s' % (i_pygccxml, prefix) )

    #execute('Unpacking Py++...', 'cd %s && unzip  %s/%s.zip' % (working_dir, data_dir, i_pyplusplus) )
    #execute('Installing Py++...', cd_work_env + 'cd %s && python setup.py install --prefix=%s' % (i_pyplusplus, prefix) )

    execute('Downloading Boost...', 'cd %s && %s "http://downloads.sourceforge.net/project/boost/boost/%s.%s.%s/%s.tar.bz2"' % ( (source_dir, Options.loader) + i_boost_version + (i_boost,) ) )
    execute('Unpacking Boost...', cd_work_env + 'tar -vjxf %s/%s.tar.bz2' % (source_dir, i_boost) )
    execute('Installing Boost...', cd_work_env + 'cd %s && ./bootstrap.sh --prefix=%s --with-libraries=python && ./bjam install --prefix=%s' % (i_boost, prefix, prefix) )

    print 'Copying Boost system library source files...',
    shutil.copytree(working_dir + '/' + i_boost + '/libs/system/src', prefix + '/include/libs/system/src')
    print 'Done!'

    #execute('Installing Boost...', cd_work_env + 'cd %s && ./bootstrap.sh --prefix=%s && ./bjam install --prefix=%s' % (i_boost, prefix, prefix) )
    # ^^^^^ alternative way to install Boost (with all libraries etc) we probably don't need that...

    print '\nSetting PyRosetta Build environment is Done!'

    print 'Creating PyRosetta build script at %s...' % i_BuildPyrosetta
    file(i_BuildPyrosetta, 'w').write(BashFileTemplate % dict(prefix=prefix, python_prefix=sys.prefix, i_python_lib=i_python_lib, compiler=Options.compiler) )
    os.chmod(i_BuildPyrosetta, stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR)
    print 'PyRosetta build script created as %s. Copy this it inside your rosetta_source dir and it run to build PyRosetta!' % i_BuildPyrosetta

    if Options.rosetta_source:
        Options.rosetta_source = os.path.abspath(Options.rosetta_source)
        print 'rosetta_source location specifying as %s. Preparing to build PyRosetta...' % Options.rosetta_source

        if not os.path.isdir(Options.rosetta_source):
            print 'Local copy of Rosetta source is not found at %s..., checking it out...' % Options.rosetta_source
            execute('Checking out fresh copy of Main rosetta/rosetta_source trunk...', 'svn co https://svn.rosettacommons.org/source/trunk/rosetta/rosetta_source %s' % Options.rosetta_source)

        print 'Copying PyRosetta build script to %s...' % Options.rosetta_source
        shutil.copy(i_BuildPyrosetta, Options.rosetta_source+'/BuildPyRosetta.sh')
        execute('Copying boost lib...', 'cp %s/lib/libboost_python.[!a]* %s/src/python/bindings/' % (prefix, Options.rosetta_source) )
        execute('Compiling PyRosetta...', 'cd %s && ./BuildPyRosetta.sh -j%s' % (rosetta_source_dir, Options.jobs))


    # we no longer need this...
    '''
    if Platform == "macos":
        print 'Platform is MacOS... adjusting gccxml config...'
        execute('Saving old gccxml config...', 'mv %s/share/gccxml-0.9/gccxml_config %s/share/gccxml-0.9/gccxml_config.original' % (prefix, prefix) )

        execute('Creating new config...', """echo 'GCCXML_COMPILER="/usr/bin/c++-4.0"' >%s/share/gccxml-0.9/gccxml_config""" % prefix)
        execute('Creating new config.....', """echo 'GCCXML_CXXFLAGS=""' >>%s/share/gccxml-0.9/gccxml_config""" % prefix)
        '''


BashFileTemplate = '''#!/bin/bash

prefix=%(prefix)s

source $prefix/PyRosettaBuildEnvironment.sh

cd src/python/bindings
./BuildBindings.py \\
    -I$prefix/include \\
    -I$prefix/include/boost \\
    -I%(python_prefix)s/include/%(i_python_lib)s \\
    --python-lib %(i_python_lib)s \\
    --boost-lib=boost_python \\
     -L . -L ./../../../../ -L$prefix/lib \\
     --compiler=%(compiler)s \\
     $*
'''


GCC_DowngradeMessage = '''Your system GCC version was reported as:\n\n%s\n\n
We currently does not support GCC 4.8* and higher. Please downgrade it!

Note: for downgrading GCC on Ubuntu Linux based system run 'sudo scripts/switch-gcc-version.sh'
'''

if __name__ == "__main__": main(sys.argv)
