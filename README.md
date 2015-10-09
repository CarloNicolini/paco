# PACO "Partitioning Cost Optimization"

Paco is a tool to build solvers and optimize quality functions on graphs for community detection.
Paco comes as both a C++ library and a support for Matlab/Octave mex function.

## Requirements
Paco needs the `C` libraries`igraph`

To compile, clone this dataset with git:

    $> git clone https:/github.com/CarloNicolini/paco.git
    $> cd paco
    $> mkdir build
    $> cd build
    $> cmake ..

If you want to compile `paco_mx` mex function, you need to enable the `MATLAB_SUPPORT` flag. At the last step:

    $> cmake -DMATLAB_SUPPORT=True ..

or equivalently for Octave

    $> cmake -DOCTAVE_SUPPORT=True ..

and then compile the source with `make`:

    $> make

## Obtaining igraph in OSX
We strongly suggest to obtain `igraph` in OSX using **Homebrew**.

    $> ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
    $> brew tap homebrew/science
    $> brew install igraph

## Obtaining igraph in Ubuntu 12.04 or greater
We suggest compiling `igraph` from source, since PACO needs `igraph`>=0.7.1.

Go to www.igraph.org/c and follow the documentation for the compilation using the autogen script they provide.
In the rest of the notes we assume that the Igraph include directory is `/usr/local/include/igraph`.

# FAQ
---

## Linking libstdc++.so.6 problem

I can compile Paco for MATLAB but after calling `paco_mx`, MATLAB prompts me with the following error message:

    Invalid MEX-file '~/paco_mx.mexa64':
    /usr/local/MATLAB/R2015a/bin/glnxa64/../../sys/os/glnxa64/libstdc++.so.6:
    version `GLIBCXX_3.4.21' not found (required by ...

This problem means that the `libstdc++.so.6` inside the Matlab library folder is pointing to a version of `libstdc++` older than the system one, usually stored in `/usr/lib/x86_64` folder.

To be sure, try:

    $> ls -al /usr/lib/x86_64/libstd*

You should get something like:

    lrwxrwxrwx 1 root root      19 Apr 23 19:00 /usr/lib/x86_64-linux-gnu/libstdc++.so.6 -> libstdc++.so.6.0.21
    -rw-r--r-- 1 root root 1541600 Apr 23 19:23 /usr/lib/x86_64-linux-gnu/libstdc++.so.6.0.21

This means that a symbolic link named `libstdc++.so.6` exists and points to the object file `libstdc++.so.6.0.21`, that is, `GLIBCXX_3.4.21`

The **solution** is to redirect the Matlab symbolic link to the system `libstdc++`.
Before doing that you'd better make a backup.
Assuming that your Matlab root folder is:

    /usr/local/MATLAB/R2015a/

then:

    $> cd /usr/local/MATLAB/R2015a/sys/os/glnxa64
    $> ls libstdc++*

and look for the object file: `libstdc++.so.6.0.XXX` (where XXX is a number like 17 in my case, since I have Matlab R2015a but can be different for previous matlab versions, this is just to specify not to look to the symbolic link), then make a backup copy:

    $> sudo cp libstdc++.so.6.0.17 libstdc++.so.6.0.17_BACKUP

and redirect the `libstdc++.so.6` to link the system library:

    $> sudo ln -fs /usr/lib/x86_64-linux-gnu/libstdc++.so.6.0.17 libstdc++.so.6

This command makes the `libstdc++.so.6` point to the system libstdc++ library.


For other informations take a look at http://dovgalecs.com/blog/matlab-glibcxx_3-4-11-not-found/
