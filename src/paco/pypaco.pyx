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
from ctypes import c_double

# Cython imports
from libcpp.string cimport string
from libcpp.map cimport map
from libcpp.vector cimport vector
import cython

ctypedef map[string, int] params_map

cdef extern from "Graph.h":
    cdef cppclass GraphC:
        GraphC() except +
        GraphC(double *A, int n, int m) except +
        void init(const double *ewlist, int _is_weighted, int m) except +
        void info() except+

cdef extern from "Community.h":
    cdef enum QualityType:
        pyQualityType
    cdef enum OptimizerType:
        pyOptimizerType
    cdef cppclass CommunityStructure:
        CommunityStructure(const GraphC *) except +
        void set_random_seed(int n)
        double optimize(QualityType quality, OptimizerType method, int repetitions)  except +
        void reindex_membership()
        vector[int] get_membership_vector()

@cython.boundscheck(False)
@cython.wraparound(False)
def paco(np.ndarray[double, ndim=2, mode="c"] graph_rep not None, **kwargs):
    """
    PACO: PArtitioning Cost Optimization
    
    Example: passing graph as adjacency matrix
      import numpy as np
      from pypaco import paco
      G = nx.karate_club_graph()
      A  = nx.to_numpy_matrix(G)
      [membership,quality] = paco(A, quality=0, nreps=10)
    
    Example: passing graph as edges list
      import numpy as np
      from pypaco import paco
      G = nx.karate_club_graph()
      E = np.array(G.edges()).astype(float64) # E must be a numpy array of doubles, not a list
      [membership,quality] = paco(E, quality=0, nreps=10)

    Example: passing graph as weighted edges list
      import numpy as np
      from pypaco import paco
      G = nx.karate_club_graph()
      E = np.array(G.edges())
      # Generate weights on the edges, in values between 0 and 10-1
      W = np.random.randint(10,size=[G.number_of_edges(),1])
      EW = np.concatenate((E,W),axis=1).astype(float)
      [membership,quality] = paco(EW, quality=2, nreps=10)

    Usage:
        [membership, quality] = paco(A, **kwargs)

    Args: 
        graph_rep: Adjacency matrix of the graph as a numpy 2D array, symmetric, binary or weighted. Otherwise an edgelist with m rows and 2 or 3 columns representing nodes indices
    Kwargs:
        quality: quality function to maximize
            0: Surprise
            1: Significance
            2: Asymptotic Surprise
            3: Infomap

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
        quality: the partition quality value
    """
    args = ['nreps','quality', 'seed', 'opt_method']

    args_diff = set(kwargs.keys()) - set(args)
    if args_diff:
        raise Exception("Invalid args:" + str(tuple(args_diff)) + "as graph_rep: valid arguments are " + str(args))

    cdef params_map par
    par[str("nreps")] = kwargs.get("nreps",1)
    par[str("quality")] = kwargs.get("quality",1)
    par[str("seed")] = kwargs.get("seed", -1)
    par[str("opt_method")] = kwargs.get("opt_method", 0)

    # Create graph instance
    cdef GraphC *G

    num_rows = graph_rep.shape[0]
    num_cols = graph_rep.shape[1]
    cdef np.ndarray[double, ndim=1, mode="c"] edges_list
    cdef np.ndarray[double, ndim=1, mode="c"] weights_list
    if num_cols<=3 and num_rows>2:
        # The matrix is an edgelist representation, passing a mx3 (if weighted) or mx2 (if binary) vector
        try:
            G = new GraphC()
            G.init(&graph_rep[0,0], int(num_cols!=2), num_rows)
        except RuntimeError:
            raise
    elif num_cols==num_rows and num_cols>2:
        # The matrix is an adjacency matrix
        try:
            G = new GraphC(&graph_rep[0,0], num_cols, num_cols)
        except RuntimeError:
            raise
    # Create community structure instance
    cdef CommunityStructure* c
    try:
        c = new CommunityStructure(G);
    except RuntimeError:
        raise 

    c.set_random_seed(int(par["seed"]));

    try:
        finalquality = c.optimize(int(par["quality"]),int(par["opt_method"]),int(par["nreps"]))
    except RuntimeError:
        raise 
    
    c.reindex_membership()
    membership = c.get_membership_vector()

    del G
    del c

    return membership, finalquality

