#   This file is part of PACO-PArtitioning Clustering Optimization a program
#   to find network partitions using modular solvers and quality functions.
#   Copyright (C) 2015 Carlo Nicolini <carlo.nicolini@iit.it>
#   PACO is free software; you can redistribute it and/or
#   modify it under the terms of the GNU Lesser General Public
#   License as published by the Free Software Foundation; either
#   version 3 of the License, or (at your option) any later version.
#
#   Alternatively, you can redistribute it and/or
#   modify it under the terms of the GNU General Public License as
#   published by the Free Software Foundation; either version 2 of
#   the License, or (at your option) any later version.
#
#   PACO is distributed in the hope that it will be useful, but WITHOUT ANY
#   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
#   FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License or the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU Lesser General Public
#   License and a copy of the GNU General Public License along with
#   PACO. If not, see <http://www.gnu.org/licenses/>.


from __future__ import division
import numpy as np
cimport numpy as np

# Cython imports
from libcpp.string cimport string
from libcpp.map cimport map
from libcpp.vector cimport vector
import cython

ctypedef map[string, int] params_map

cdef extern from "Graph.h":
    cdef cppclass GraphC:
        GraphC(double *A, int n, int m)

cdef extern from "Community.h":
    cdef enum QualityType:
        pyQualityType
    cdef enum OptimizerType:
        pyOptimizerType
    cdef cppclass CommunityStructure:
        CommunityStructure(const GraphC *)
        void set_random_seed(int n)
        double optimize(QualityType quality, OptimizerType method, int repetitions)
        void reindex_membership()
        vector[int] get_membership_vector()

@cython.boundscheck(False)
@cython.wraparound(False)
def paco(np.ndarray[double, ndim=2, mode="c"] input not None, **kwargs):
    """
    PACO: PArtitioning Cost Optimization
    
    Example:
      import numpy as np
      from pypaco import paco
      G = nx.karate_club_graph()
      A  = nx.to_numpy_matrix(G)
      [membership,quality] = paco(A, quality=0, nreps=10)
    
    Usage:
        [membership, quality] = paco(A, **kwargs)

    Args: 
        input: Adjacency matrix of the graph as a numpy 2D array, symmetric, binary or weighted.
    Kwargs:
        quality: quality function to maximize
            0: Surprise
            1: Significance
            2: Asymptotic Surprise
            3: Infomap
            4: Modularity (EXPERIMENTAL)
            5: Asymptotic Modularity (EXPERIMENTAL)

        opt_method: Optimization method
            0: Agglomerative,
            1: Random,
            2: SimulatedAnnealing,
            3: Infomap

        nreps: number of repetitions of PACO (can increase the quality of the partition), (default 1).
        seed: random seed for randomization (integer value)
        
        save_solution: 1 to save ongoing solution, 0 otherwise (default 0)
        
    Out:
        membership: a list of vertices community membership
        surprise: the value of Surprise for the current partition
    """
    args = ['nreps','quality', 'seed', ',opt_method']

    args_diff = set(kwargs.keys()) - set(args)
    if args_diff:
        raise Exception("Invalid args:" + str(tuple(args_diff)) + "as input: valid arguments are " + str(args))

    cdef params_map par
    par[str("nreps")] = kwargs.get("nreps",1)
    par[str("quality")] = kwargs.get("quality",1)
    par[str("seed")] = kwargs.get("seed", -1)
    par[str("opt_method")] = kwargs.get("opt_method", 0)

    n = input.shape[0]
    cdef GraphC *G = new GraphC(&input[0,0], n, n)
    cdef CommunityStructure* c = new CommunityStructure(G);
    c.set_random_seed(int(par["seed"]));

    finalquality = c.optimize(int(par["quality"]),int(par["opt_method"]),int(par["nreps"]))
    c.reindex_membership()
    membership = c.get_membership_vector()

    del G
    del c

    return membership, finalquality

