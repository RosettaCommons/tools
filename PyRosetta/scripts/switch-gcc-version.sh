#!/bin/bash

apt-get install gcc-4.1 g++-4.1 python-all-dev libzip-dev subversion cvs

update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.1  310
update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.2  210
update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.3  100
update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.4  180

update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-4.1  310
update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-4.2  210
update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-4.3  100
update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-4.4  180


update-alternatives --install /usr/bin/cc  cc /usr/bin/gcc 100
update-alternatives --set cc /usr/bin/gcc

update-alternatives --install /usr/bin/c++  c++ /usr/bin/g++ 100
update-alternatives --set c++ /usr/bin/g++


echo "Settings done, now configure using autoconfig..."
update-alternatives --config gcc
update-alternatives --config g++
