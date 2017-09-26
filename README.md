# Experimental GR Particle Pusher

## How to compile

Checkout the repository with `git`:

    git clone https://github.com/fizban007/gr_pusher.git

Go into the new `gr_pusher` directory. There, create a new subdirectory called `build`:

    mkdir build
    cd build
    cmake ..

All external dependencies are bundled in the repo, so as long as you have a
decently recent compiler that has `c++14` support, it should compile. This means
GCC 5.0+, Clang 3.4+, or Intel C++ compiler 17+. The code has been tested mainly
on Intel parallel studio 2017 and GCC 7.1.1.

To compile the main program, simply run

    make

in the `build` directory. To run the (very incomplete) test-suite, run

    make check

in the same directory. Note that GCC for whatever reason will take significantly
longer than Intel compiler to compile the entire test suite, due to heavy use of
template metaprogramming.

After compiling, some executable starting with `pusher` should be in the
`gr_pusher/bin` directory. Run it to execute the test case that is written in
`main.cpp`.
