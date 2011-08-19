This python script uses clang to generate an XML version of the ASTs
and then determine which functions access core::options::option[]. To
use these functions, you have to build clang yourself, which isn't that
difficult. Instructions for installing it are given here:

http://clang.llvm.org/get_started.html

The instructions are slightly out of date and create a Debug+Assert build,
which is slow and throws an assertion when dumping the XML file. To
create a release build without assertions, do the following:

svn co http://llvm.org/svn/llvm-project/llvm/trunk llvm
cd llvm/tools
svn co http://llvm.org/svn/llvm-project/cfe/trunk clang
cd ..
./configure --enable-optimized --disable-assertions --with-cxx-include-root=/usr/include/c++/4.4.0 --with-cxx-include-arch=x86_64-redhat-linux6E --with-cxx-include-32bit-dir=32
make -j 5

Note that clang requires libstdc++ version 4.2 or higher. I used the
gcc44-c++ and libstdc++44-devel packages from CentOS to install them
on my linux workstation. You need to tell clang where the headers are,
as shown above.

Before running the script, the directory with the clang binary must
be in your PATH environment variable. You may also need to change the
platform path if you're not using linux.

The script expects to be run from the src directory. After configuring
is as described above, simply run the script without arguments. If you
want to have it print the uses of options as it finds them, and not just
at the end, uncomment the last print statement in the script.
