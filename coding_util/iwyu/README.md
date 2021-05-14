Include What You Use For Rosetta
================================

Include What You Use (IWYU: <https://include-what-you-use.org/>) is a Clang-based tooling
which attempts to optimize which headers are included in your C++ files.
(That is, you include only the headers needed for compilation, but you include all such headers, even transitive ones.)

The default IWYU approach doesn't really work as-is with the way Rosetta approaches things,
so the scripts in this directory are intended to manipulate the IWYU output in order to make it more "Rosetta".

Like the standard IWYU approach, the scripts are divided into two phases, run and apply.

IWYU Version
------------

IWYU is currently under heavy development, and the fixes it suggests depend in part on version.

The scripts in this directory are currently built to work with version 0.14 (Clang version 10.0.0).
They may or may not work properly with newer/older versions of IWYU, 
and some updates to the approach may be needed for newer IWYU versions.

Run
---

The `run_iwyu.py` script handles the run phase. 
The intent is that `run_iwyu.py` script will be run against a clean (no working directory changes) repo,
and will generate `*.riwyuf` files ('Rosetta include what you use fixes') files describing
an edited collection of the suggested fixes.

*Suggested Run Command*

The repository must be in a compilable state prior to running IWYU over it.

While you can run on individual files, it's recommended to run on the entire repo
before applying any edits. (As edits to one file can result in other files then not compiling.)

    cd Rosetta/main/source/
    ../tools/coding_util/iwyu/run_iwyu.py  src/ test/

There may be some error and warning messages. 
It's best to fix those issues and re-run before moving on to the apply stage.

*Differences from  default*

The main job of the script is to run IWYU on each hh/cc file, and parse the results.
There's several edits made to the results. 
The biggest of which is to convert explicit forward declarations into inclusion of forward headers.
There's also some alterations to account for other foibles of differences in how IWYU and Rosetta handle things.

*Support files*

`Rosetta.imp` -- IWYU has an input file which describes alterations which should be applied.
Rosetta.imp is this IWYU-formatted (see the official documentation) input file fed directly to the IWYU applicaiton.
As much as possible, fixes for Rosetta should be represented in that fashion 

`IWYU_nonstandard_fwd.txt` -- Most of the conversions of explicit forward declarations suggested by IWYU
into forward header inclusion are straightforward.
But some forward declarations are in header names which are not autodiscoverable. The `IWYU_nonstandard_fwd.txt` lists these. 
Needed additions to this file should be obvious, as the `run_iwyu.py` script should print
`WARNING: Can't find forward header file:` for potnetial files.
NOTE: Please don't list regular hh files here - it's almost always worth making a fwd.hh file if one doesn't already exists.

`IWYU_forced_subs.txt` -- IWYU sometimes wants to expose internal implementation files. 
The forced sub file lists internal implementation files and the public file which should be included instead.
(Note that the Rosetta.imp should theoretically take care of this, but it doesn't always.)

*Output*

The `*.riwyuf` files are JSON formatted files with the edited suggested fixes, 
and will be written in the same directory as the coresponding source file.

Apply
-----

A naive application of all of the fixes listed in the files (even with the Rosetta-specific edits) doesn't always work,
due to 1) bugs in IWYU and 2) differences in how Rosetta handles things versus how IWYU thinks things should be handled.
As such, the `apply_iwyu.py` script takes the approach of test compilation. 
That is, it applies fixes one-by-one, and checks for compilation of the file.

The approach `apply_iwyu.py` takes is to prioritize removal/deletions of header files, 
and add in the minimal number of headers to make sure the file compiles. 
This is slightly different than the IWYU approach, in that it won't necessarily add transitive headers.
As such, you may need to run `apply_iwyu.py` multiple times in order to correct fully for header removals.

