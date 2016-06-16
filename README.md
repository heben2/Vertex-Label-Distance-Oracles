To run on Linux (tested on Arch)
Install
- Cmake
- c++ x14
Run commands
    cd src/
    cmake .
    make
    ./main

The Makefile is included, so might work out of the box with make from /src/ on 
unix systems with make.

To run on Windows (tested on win7)
Install 
- Visual Studio 14 2015
- CMake
Install with cmake (terminal/cmd):
    cd src/
    mkdir build     (can be done manually, i.e. make a new directory called 'build')
    cd build
Then run command
    cmake -G "Visual Studio 14 2015" ..

Now import the project in Visual Studio 14 2015, compile (F7) 
run src/build/Release/main.exe


The oracles are located in the oracle library is located in /oracles/
The tests methods are located in /tests/
Note that the tests methods are quite a mess and not made to be run without 
editing.
