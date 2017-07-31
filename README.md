Rosetta Tools
=============

Rosetta tools contains a number of utilities to facilitate Rosetta development including tools for C++ parsing, header only compilation and measuring test coverage.  Somewhat mysteriously, this repository also contains scripts to run particular Rosetta applications including ERRASER and fragment picking.

Protein Tools
=============

A Collection of useful python modules that are used by many scripts here. 

## Installation:


The library depends on Biopython
To use the graphics package, you will also need to install matplotlib

```
pip install biopython --user
pip install matplotlib --user
```

The library is divided into 3 packages:

graphics - Classes and functions for graphical output (requires matplotlib)
rosetta - Classes and functions for parsing file formats specific to rosetta
protein - Classes and functions for manipulation protein structural data (requires biopython)

to install this package type:


### Typical Install
```
cd protein_tools
python setup.py install
```

### Mac Install
```
cd protein_tools
python setup.py install --user
```

### Custom Install
```
cd protein_tools
python setup.py install --install-scripts=/script/install/directory
```
