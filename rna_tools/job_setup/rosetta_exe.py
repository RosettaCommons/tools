from os.path import exists, join
import os
from itertools import product

def rosetta_exe(exe_file, rosetta_folder = None, extension = None) :
    """
    Return the absolute path to the given exexutable file in Rosetta.
    If no input is given, return the Rosetta bin folder path.
    Only the prefix should be given. Rosetta will figure out the rest.
    User can provides either the Rosetta root path or the binary path as optional input,
    or the function will search for $ROSETTA environment variable.
    """
    if rosetta_folder is None:
        try:
            rosetta_folder = os.environ["ROSETTA"]
        except KeyError:
            pass
    name_extensions = ['', ".linuxgccrelease", ".linuxclangrelease", ".macosgccrelease", ".macosclangrelease"]

    if rosetta_folder is not None:
        exe_folder = join(rosetta_folder, "main/source/bin/") #Default Rosetta folder structure
        if not exists(exe_folder): #Otherwise, assume the input folder name is bin path
            exe_folder = rosetta_folder
        check_path_exist(exe_folder)
        exe_path = ""
        if extension is not None:
            if extension[0] != '.':
                extension = '.'+extension
            exe_path = join(exe_folder, exe_file + extension)
        else:
            for ext in name_extensions :
                exe_path = join(exe_folder, exe_file + ext)
                if exists(exe_path) :
                    break
        check_path_exist(exe_path)
        return exe_path
    else:
        paths = os.environ['PATH'].split(':')
        if extension is not None:
            for path in paths:
                exe_path = join(path, exe_file + ext)
                if exists(exe_path):
                    break
        else:
            for path, ext in product(paths, name_extensions):
                exe_path = join(path, exe_file + ext)
                if exists(exe_path) :
                    break
    check_path_exist(exe_path)
    return exe_path

#####################################
def check_path_exist(path_name) :
    if not exists(path_name) :
        raise ValueError("Path %s does not exist!" % path_name)
