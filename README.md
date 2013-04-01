Rosetta Tools
=============

Rosetta tools contains a number of utilities to facilitate Rosetta development including tools for C++ parsing, header only compilation and measuring test coverage.  Somewhat mysteriously, this repository also contains scripts to run particular Rosetta applications including ERRASER and fragment picking.

Configure Rosetta
-----------------
via `curl`
```
curl -u $USERNAME -L https://github.com/RosettaCommons/rosetta_tools/raw/master/configure_rosetta_repo.sh -o rosetta.sh | sh rosetta.sh
```

via `wget`

```
wget --no-check-certificate https://github.com/RosettaCommons/rosetta_tools/raw/master/configure_rosetta_repo.sh -O - | sh
```

