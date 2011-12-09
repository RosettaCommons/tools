Intro
=====

Motivation
----------

I am tired of writing complicated regular expressions to try and automate C++ refactoring. There is no fundamental reason we have to spend hundreds of man-hours fixing compile errors at XRW. It’s time for some C++ parsing real tools.

This is an attempt at that based on the clang compiler. ATM, it's a bit slow compared to compiling (or grep), but with libclang, it will (eventually) be possible to do refactoring operations *completely* and *correctly* in a turn-key fashion.

.. todo::
	improve speed... removing necessity for parent links and locmap will help

What already works
------------------

You can

* completely parse C++ source files
* get the full namespace of anything, regardless of context
* get the headers that a file *really* needs


Requirements
------------

Libclang
~~~~~~~~

I got the following from directions here: http://clang.llvm.org/get_started.html#build

mkdir libclang

cd libclang/

svn co http://llvm.org/svn/llvm-project/llvm/trunk llvm

cd llvm/tools/

svn co http://llvm.org/svn/llvm-project/cfe/trunk clang

cd ../..

mkdir build

cd build/

../llvm/configure 

make -j8 ENABLE_OPTIMIZED=1

# clang is here: build/Release+Asserts/bin/clang++ 

# libclang is here: build/Release+Asserts/lib

export LD_LIBRARY_PATH=/PATH_TO_LIBCLANG/build/Release+Asserts/lib:LD_LIBRARY_PATH

export PYTHONPATH=/PATH_TO_LIBCLANG/llvm/tools/clang/bindings/python:PYTHONPATH


Sphinx
~~~~~~

for documentation, you need sphinx, which is the defacto python doc tool.

sudo pip install sphinx

or, if you don't have pip 'sudo easy_install sphinx', but you should really really get pip (sudo easy_install pip will do it)

anything checked in really should have up-to-date docs, but you can rebuild docs with:

cd doc; make html

TBD unittest framework?
~~~~~~~~~~~~~~~~~~~~~~~