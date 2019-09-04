import os
import sys
import setuptools

def setup_package():
    # chdir to source path and add sources to sys.path to allow invocation of setup script from other directories.
    # Code structure taken from pyrosetta/setup/setup.py, which got it from numpy/setup.py
    src_path = os.path.dirname(os.path.abspath(sys.argv[0]))
    old_path = os.getcwd()
    os.chdir(src_path)
    sys.path.insert(0, src_path)

    try:
        setuptools.setup(
            name = 'ambrose',
            version = '2.1',
            author = 'Maria Szegedy',
            author_email = 'mszegedy2@gmail.com',
            description = 'AMBER interoperability for PyRosetta.',
            url = 'https://github.com/RosettaCommons/tools/tree/master/AmbRose',
            packages = setuptools.find_packages(exclude=['legacy']),
            package_data = {
                '' : ['*.prototype']
            },
            classifiers = [
                'Programming Language :: Python :: 3'
            ]
        )
    finally:
        os.chdir(old_path)
        del sys.path[0]

if __name__ == '__main__':
    setup_package()
