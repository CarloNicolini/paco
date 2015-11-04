# PACO "Partitioning Cost Optimization"

PACO is a tool to build solvers and optimize quality functions on graphs for community detection. PACO comes as both a C++ library and a support for Matlab/Octave mex function.

At this moment PACO allows to optimize the following quality functions:

1. Newman's Modularity (binary and weighted networks)
2. Surprise (only binary networks)
3. Asymptotical Surprise (binary and weighted networks)
4. Significance (binary and weighted networks)

and supports three optimization algorithms:

1. An agglomerative optimizer based on FAGSO, a very fast and powerful heuristics that is based on sorting of networks edges by their endpoints neighbors set Jaccard coefficient.
2. Simulated Annealing optimizer: the classical simulated annealing algorithm for global optimization.
3. Random optimizer

## Requirements
PACO is written in C++ and therefore needs a viable C/C++ compiler, typicall `gcc` or `clang` are supported compilers. PACO is tested to compile on Linux and OSX.

PACO makes deep use of the `C` libraries`igraph`. To see how to get `igraph` look at the appropriate section .

To get and compile PACO, issue the following commands in the command-line:

    $> git clone https:/github.com/CarloNicolini/PACO.git
    $> cd PACO
    $> mkdir build
    $> cd build
    $> cmake ..

If you want to compile `PACO_mx` mex function, you need to enable the `MATLAB_SUPPORT` flag. At the last step:

    $> cmake -DMATLAB_SUPPORT=True ..

or equivalently for Octave

    $> cmake -DOCTAVE_SUPPORT=True ..

and then compile the source with `make`:

    $> make

# Lancichinetti-Fortunato-Radicchi benchmark functions
Lancichinetti-Fortunato-Radicchi benchmark is a network generator for benchmark of community detection algorithms.

In this release I've adapted and partially rewritten the freely-available code (https://sites.google.com/site/santofortunato/inthepress2) to generate weighted networks as described in the paper by Lancichinetti A. and Fortunato S., *Benchmarks for testing community detection algorithms on directed and weighted graphs with overlapping communities*, Phys. Rev. E 80, 016118.

This release of PACO comes with some convenience MATLAB wrappers around the famous LFR benchmark. You can compile the LFR benchmark by specifying:

    $> cmake -DCOMPILE_LFR=True -DMATLAB_SUPPORT=True ..
    $> make lfrw_mx

The `lfrw_mx` function is self-documented:

    LFRW: Lancichinetti-Fortunato-Radicchi network generator for weighted networks
    This program produces weighted adjacency matrices for graph community detection benchmark
    
    [WeightedGraph, membership] = lfrw_mx('argname',argvalue);
    
    Input:
    'N'
        Desidered number of nodes in the network
    'k'
        Desidered average degree of the network
    'maxk'
        Desidered max degree of the network
    'mut'
        Desidered community topological mixing coefficient (range [0,1])
    'muw'
        Desidered weights mixing coefficient (range [0,1])
    't1'
        Desidered exponent for the degree distribution
    't2'
        Desidered exponent for the community size distribution
    'beta'
        Desidered beta exponent
    'C'
        Desidered clustering coefficient

An example usage of `lfrw_mx` in Matlab is the following:

1. Start a MATLAB console (no graphical interface for simplicity):
```
$> matlab -nojvm
```

2. Once MATLAB is started:
```
>> [AdjMatrix, planted_vert_membership] = lfrw_mx('N',100,'k',5,'maxk',15,'muw',0.2,'mut',0.2);
```

3. Some informations are printed on output, regarding the generation of the benchmark network. The first output argument is the adjacency matrix of the generated network, the second is the planted partition indicated as a vector of integers, where every element is the community of the i-th vertex.

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

I can compile PACO for MATLAB but after calling `PACO_mx`, MATLAB prompts me with the following error message:

```
Invalid MEX-file '~/PACO_mx.mexa64': /usr/local/MATLAB/R2015a/bin/glnxa64/../../sys/os/glnxa64/libstdc++.so.6: version `GLIBCXX_3.4.21' not found (required by ...
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
