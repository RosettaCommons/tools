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

Rosetta Version
---------------

Both `run_iwyu.py` and `apply_iwyu.py` rely on having the current Clang++ compile line. 
This is encoded in the `iwyu_support.py` file.

Setup
-----

Include what you use should be installed into a standard location, 
ideally such that the default `iwyu` command will invoke it.
The corresponding Clang compiler should be invokable with the default `clang++`.
(Although there should be options on the scripts to change these, if necessary.)

Suggested Usage
---------------

Feel free to vary from this usage, as you preference dictates.


1. Start with a clean working directory which compiles with Clang
2. Run `run_iwyu.py` on the full src/ and test/ directories
3. Fix any issues mentioned (commiting fixes), and repeat `run_iwyu.py` until you get a clean run.
4. Run `apply_iwyu.py`
5. Keep re-running `apply_iwyu.py` until the script reports no modifications.
6. Do a `grep -R 'AUTO IWYU' src/ test/` to check for over-zealous additions
7. Check compilation by normal means, making manual fixes, as needed
    Manual fixes may be "just check out the original version of some files"
8. Once you get a clean compile, run `sed -i -r '/AUTOREMOVED IWYU/d' $(grep -Rl 'AUTOREMOVED IWYU' src/ test/)`
    to remove the autoremoved lines.
9. Double check you're still a clean compile
10. Commit changes to src/ and test/
11. (Optional) Update the `run_iwyu.py` and `apply_iwyu.py` scripts to handle failure cases encounterd.
12. Repeat steps 2-10 until the process no longer suggests changes.

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
    ../tools/coding_util/iwyu/run_iwyu.py -j12 src/ test/

Alter the -j setting for the number of processors availible.

There may be some error and warning messages. 
It's best to fix those issues and re-run before moving on to the apply stage.

Note that due to the way IWYU wants to handle forwarding, the script won't run on fwd.hh by default.
Those are best handled manually.

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
`WARNING: Can't find forward header file:` for potential files.
NOTE: Please don't list regular hh files here - it's almost always worth making a fwd.hh file if one doesn't already exists.

`IWYU_forced_subs.txt` -- IWYU sometimes wants to expose internal implementation files. 
The forced sub file lists internal implementation files and the public file which should be included instead.
(Note that the Rosetta.imp should theoretically take care of this, but it doesn't always.)

*Common Errors/Warnings*

`WARNING: Can't find forward header file:` -- Add a forward header for the given symbol. 
If a forward header for the symbol already exists, you can add the translation to `IWYU_nonstandard_fwd.txt` 

`WARNING: File <fn> does not compile with regular Clang.` -- As the script doesn't use the regular build system,
it can find files which aren't compiled normally. As such, there's a check of compilation against plain clang
if the IWYU errors out. If you don't expect these files to compile (e.g. they're non-compiled pilot or devel apps)
then there's nothing to be done. If you expect them to compile, though, you should figure out if there's compilation issues.
(Note that you'll also get this error if the compilation backcheck is not working for some reason.)

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

*Suggested Run Command*

    cd Rosetta/main/source/
    ../tools/coding_util/iwyu/apply_iwyu.py -j12 src/ test/

*Repeated application*

Since removing a header can have a downstream effect on the compile-ability of other files,
it's recommended to keep running the `apply_iwyu.py` command until it reports that no changes to files have been made.
(Running `apply_iwyu.py` can fix the compilation of files if the reason for the compilation is a missing header
which is listed in the corresponding .riwyuf file. While `run_iwyu.py` can't generate a riwyuf file from a non-compiling file,
if you have have the version from when it was compiling, it may be able to fix things.)

Note that subsequent runs of `apply_iwyu.py` after the first one are much faster. (And can/should be applied without resetting the file contents.)

Some issues may not be fixable with repeated applications of `apply_iwyu.py`.
You will likely need to try compiling and fix (hopefully a scant few) things up manually.

*Reported issues*

Unlike the `run_iwyu.py` script, issues with `apply_iwyu.py` don't necessarily need to be fixed prior to accepting the result.
The script should leave the source files in a compilable state.

If you do want to address the issues, I recommend doing it in a separate step.

*Fixing Common Issues*

One common issue is that IWYU gets overzealous, even with the compile check. (Particularly the case with the test/ directory.)
This will show up in files that don't compile after the run (even after repeated calls). 
Fix this by uncommenting the needed header, and then adding `// DO NOT AUTO-REMOVE` to the line

IWYU can also add some odd headers, or a somewhat extensive set of headers.
The scripts are set up to avoid most of that, but it sometimes happens.
I'd recommend doing a `grep -R 'AUTO IWYU'` to see what headers are being added (particularly boost and stdlib headers),
and which files have a large number of low-level files being added.
Often just checking out the original version and re-running will fix the odd additions.

`ERROR - Line to remove does not match expected content` -- This generally happens on subsequent rounds 
when the auto added headers are inserted above headers to remove. The inserting at the top is typically due to having non-include lines
interspersed within the include lines. (It's recommended to fix those for subsequent iterations of run/apply.)

Namespaces: IWYU 0.14 has a weird way of ascribing namespaces, so `run_iwyu.py` tries to ignore for-namespace-only includes.
This means that compilation can get broken due to `using namespace` directives which aren't actually used.
The fix is to manually delete the useless using directives.

Common reasons for "Issues Removing"
* Extra "using namespace" directives.

