Rosetta Tools
=============

Rosetta tools contains a number of utilities to facilitate Rosetta development including tools for C++ parsing, header only compilation and measuring test coverage.  Somewhat mysteriously, this repository also contains scripts to run particular Rosetta applications including ERRASER and fragment picking.

Configure Rosetta
-----------------
via `curl`: Be sure to replace `$USER` with your GitHub user name!
```
tmp=$(date +%Y%m%d%H%M); curl -u $USER -L https://github.com/RosettaCommons/rosetta_tools/raw/master/configure_rosetta_repo.sh > $tmp && sh $tmp && rm $tmp
```
