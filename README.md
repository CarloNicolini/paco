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

Go to www.igraph.org/c and follow the documentation for the compilation using the autogen script they provide, otherwise follow the following instructions:

### Compilation and installation of igraph on Ubuntu 12.04 or greater

1. Download igraph-0.7.1 source code:

```
$> cd
$> wget http://igraph.org/nightly/get/c/igraph-0.7.1.tar.gz
```
 
2. Extract the code under your home folder:

```
$> tar -zxvf igraph-0.7.1.tar.gz
```

3. Install the necessary libs, under Ubuntu:

```
$> sudo apt-get install build-essential libxml2-dev zlib1g-dev
```

4. Try to configure and compile igraph:

```
$> cd igraph-0.7.1/
$> ./configure
$> make 
```

5. If you have many processors on your computer, make in parallel mode, instead of `make` just issue `make -j8` for example to compile with 8 cores.

6. Check if installation went fine:

```
$> make check
```

7. If tests are passed, then install it. I usually install *igraph* under `/usr/lib/`. To do that I call `make install` with administrator privileges. 

```
$> sudo make install
```

8. If you don't have administrator privileges you can compile igraph under your home folder.

In the rest of the notes we assume that the Igraph include directory is `/usr/local/include/igraph`.


# Windows support:
Despite everything should be ready to be ported happily in Windows, I don't have time to let the code compile smoothly on Windows. If you want to join me in writing PACO for Windows let me know.

# FAQ
---

## Linking libstdc++.so.6 problem

I can compile Paco for MATLAB but after calling `paco_mx`, MATLAB prompts me with the following error message:

```
Invalid MEX-file '~/paco_mx.mexa64': /usr/local/MATLAB/R2015a/bin/glnxa64/../../sys/os/glnxa64/libstdc++.so.6: version `GLIBCXX_3.4.21' not found (required by ...
```

This problem means that the `libstdc++.so.6` inside the Matlab library folder is pointing to a version of `libstdc++` older than the system one, usually stored in `/usr/lib/x86_64` folder.

To solve the issue you need to redirect the symbolic links in the MATLAB folder to the systemwise `libstdc++`. Hereafter we assume the MATLAB folder to be `/usr/local/MATLAB/R2015a` and the system to be some Linux variant.

Two of the symlinks for libraries need to be changed:

```
$> cd /usr/local/MATLAB/R2015a/sys/os/glnxa64
$> ls -l
```

The sym links for libstdc++.so.6 and libgfortran.so.3 should point to versions in /usr/lib, not local ones.


Before changing this libraries, first make sure `g++-4.4` and `libgfortran3`are installed :

```
$> sudo apt-get install g++-4.4 libgfortran3
```

Now, modify the symlinks:

```
$> sudo ln -fs /usr/lib/x86_64-linux-gnu/libgfortran.so.3.0.0 libgfortran.so.3
$> sudo ln -fs /usr/lib/gcc/x86_64-linux-gnu/4.4/libstdc++.so libstdc++.so.6
```

This command makes the `libstdc++.so.6` point to the `g++-4.4` `libstdc++` library.

For other informations take a look at http://dovgalecs.com/blog/matlab-glibcxx_3-4-11-not-found/ or https://github.com/RobotLocomotion/drake/issues/960 which are very similar problems.
