/* This file is part of PACO-PArtitioning Clustering Optimization a program
* to find network partitions using modular solvers and quality functions.
*
*  Copyright (C) 2016 Carlo Nicolini <carlo.nicolini@iit.it>
*
*  PACO is free software; you can redistribute it and/or
*  modify it under the terms of the GNU Lesser General Public
*  License as published by the Free Software Foundation; either
*  version 3 of the License, or (at your option) any later version.
*
*  Alternatively, you can redistribute it and/or
*  modify it under the terms of the GNU General Public License as
*  published by the Free Software Foundation; either version 2 of
*  the License, or (at your option) any later version.
*
*  PACO is distributed in the hope that it will be useful, but WITHOUT ANY
*  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
*  FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License or the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU Lesser General Public
*  License and a copy of the GNU General Public License along with
*  PACO. If not, see <http://www.gnu.org/licenses/>.
*
* If you use PACO for you publication please cite:
*
* "Modular structure of brain functional networks: breaking the resolution limit by Surprise"
* C. Nicolini and A. Bifone, Scientific Reports
* doi:10.1038/srep19250
*
* "Community detection in weighted brain connectivity networks beyond the resolution limit", 
* C.Nicolini. C.Bordier, A.Bifone, arxiv 1609.04316
* https://arxiv.org/abs/1609.04316
*/


#include <iostream>
#include <string.h>
#include <sstream>
#include <igraph.h>
#ifdef UNIX
#include <sys/time.h>
<<<<<<< HEAD
#endif
#ifdef WIN32
	#ifndef strcasecmp
		#define strcasecmp _stricmp
	#endif
#endif
=======
#include <unistd.h>

>>>>>>> bbc4bcf27a1f1420de2a1129702bf4d8aff4bd97
#include "igraph_utils.h"
#include "Graph.h"
#include "Community.h"

#include "RandomOptimizer.h"
#include "AgglomerativeOptimizer.h"
#include "AnnealOptimizer.h"

#include "SurpriseFunction.h"
#include "AsymptoticSurpriseFunction.h"
#include "SignificanceFunction.h"

#ifdef __linux__
#include <mex.h>
#endif

#ifdef __APPLE__
#include "mex.h"
#endif

<<<<<<< HEAD
#ifdef WIN32
#include <mex.h>
#endif
=======
#include "mexInterrupt.h"
>>>>>>> bbc4bcf27a1f1420de2a1129702bf4d8aff4bd97

typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> MatrixXdCol;
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> MatrixXdRow;
using namespace std;

void printUsage()
{
    mexPrintf("PACO PArtitioning Cost Optimization.\n");
    mexPrintf("[membership, qual] = paco(W);\n");
    mexPrintf("Input:\n");
    mexPrintf("	W: an undirected weighted network with positive edge weights. Negative edge weights are not seen as edges and therefore discarded. Remember to use real matrices, logical matrices throw error.\n");
    mexPrintf("Output:\n");
    mexPrintf("	membership: the membership vector that represents the community which every vertex belongs to after quality optimization.\n");
    mexPrintf("	qual: the current quality of the partition.\n");
    mexPrintf("\n");
    mexPrintf("Options:\n");
    mexPrintf("paco accepts additional arguments to control the optimization process\n");
    mexPrintf("[m, qual] = paco(W,'method',val);\n");
    mexPrintf("	val is one of the following integers: {0,1,2}:\n");
    mexPrintf("		0: Agglomerative\n");
    mexPrintf("		1: Random\n");
    mexPrintf("		2: Annealing (EXPERIMENTAL)\n");
    mexPrintf("[m, qual] = paco(W,'quality',val);\n");
    mexPrintf("	val is one of the following integers: {0,1,2,3}:\n");
    mexPrintf("		0: Surprise (discrete)\n");
    mexPrintf("		1: Significance\n");
    mexPrintf("		2: AsymptoticSurprise\n");
    mexPrintf("		3: Infomap\n");
#ifdef EXPERIMENTAL_FEATURES
    mexPrintf("		4: Modularity\n");
    mexPrintf("		5: AsymptoticModularity (EXPERIMENTAL)\n");
    mexPrintf("		6: Wonder (EXPERIMENTAL)\n");
#endif
    mexPrintf("[m, qual] = paco(W,'nrep',val)\n");
    mexPrintf("	val is the number of repetitions to run over which to choose the best quality value (the lowest for Infomap, the highest for the other methods\n");
    mexPrintf("[m, qual] = paco(W,'seed',val)\n");
    mexPrintf(" val is a specific random seed to the algorithm, in order to have reproducible results.\n");
    mexPrintf("\n\n");
    mexPrintf("Example:\n");
    mexPrintf("%Create a random symmetric thresholded network\n");
    mexPrintf(">> A=rand(100,100); A=(A+A')/2; A=A.*(A>0.5);\n");
    mexPrintf("\t % Run Asymptotical Surprise optimization on A for 1000 repetitions and return the highest Surprise\n partition membership together with the quality value\n");
    mexPrintf(">> [memb, qual] = paco(A,'method',2,'nrep',1000);\n");
}

enum error_type
{
    PACO_NO_ERROR = 0,
    ERROR_TOO_MANY_OUTPUT_ARGS = 1,
    ERROR_NOT_ENOUGH_ARGS = 2,
    ERROR_ARG_VALUE = 3,
    ERROR_ARG_TYPE = 4,
    ERROR_MATRIX = 5,
    ERROR_ARG_EMPTY=6,
    ERROR_ARG_UNKNOWN=7
};

static const char *error_strings[] =
{
    "",
    "Too many output arguments.",
    "Not enough input arguments.",
    "Non valid argument value.",
    "Non valid argument type.",
    "Non valid input adjacency matrix. PACO accepts symmetric real dense-type (n x n) matrices or sparse edges-list representation \
    [num_edges x 3] array of edges list with edge endpoints and weight.",
    "Expected some argument value but empty found.",
    "Unkwown argument."
};

struct PacoParams
{
    OptimizerType method;
    QualityType qual;
    size_t nrep;      // Maximum number of consecutive repetitions to perform.
    int rand_seed; // random seed for the louvain algorithm
    int verbosity_level;
};

error_type parse_args(int nOutputArgs, mxArray *outputArgs[], int nInputArgs, const mxArray * inputArgs[], PacoParams *pars, int *argposerr )
{
    if (nInputArgs < 1)
    {
        *argposerr = 0;
        return ERROR_NOT_ENOUGH_ARGS;
    }

    if (nOutputArgs>2)
    {
        *argposerr = 0;
        return ERROR_TOO_MANY_OUTPUT_ARGS;
    }

    const mxArray *W = inputArgs[0];

    int M = mxGetM(W);
    int N = mxGetN(W);
    // In this case we are feeding instead of the adjacency matrix, the result of [i j w]=find(A);
    bool feedingSparseMatrix=false;
    if (N==3 && M>3)
    {
        feedingSparseMatrix=true;
    }

    bool v1 = M!=N;
    if (feedingSparseMatrix)
        v1=false;
    bool v2 = mxIsComplex(W);
    bool v3 = mxIsEmpty(W);
    bool v4 = mxIsCell(W);
    bool v5 = !mxIsNumeric(W);
    //bool v6 = !mxIsSparse(W);

    if ( v1 || v2 || v3 || v4 || v5 )
    {
        *argposerr = 0;
        return ERROR_MATRIX;
    }

    // Iterate on function arguments
    int argcount=1;
    while (argcount<nInputArgs)
    {
        // Be sure that something exists after c-th argument
        if (argcount+1 >= nInputArgs)
        {
            *argposerr = argcount;
            return ERROR_ARG_EMPTY;
        }
        // Couple argument type - argument value
        const mxArray *partype = inputArgs[argcount];
        const mxArray *parval = inputArgs[argcount+1];
        char* cpartype;
        // To be a valid parameter specification it must be a couple ['char',real]
        if (mxIsChar(partype) && !mxIsChar(parval))
        {
            cpartype = mxArrayToString(partype);
            //mexPrintf("ARGUMENT: %s VALUE=%g\n", cpartype,*mxGetPr(parval));
            // Parse string value inputArgs[c]
            if ( strcasecmp(cpartype,"Method")==0 )
            {
                pars->method = static_cast<OptimizerType>((int)*mxGetPr(parval));
                if (pars->method<0 || pars->method>3)
                {
                    *argposerr = argcount+1;
                    return ERROR_ARG_VALUE;
                }
                argcount+=2;
            }
            else if ( strcasecmp(cpartype,"Quality")==0 )
            {
                pars->qual = static_cast<QualityType>((int)*mxGetPr(parval));
                if (pars->qual<0 || pars->qual>3)
                {
                    *argposerr = argcount+1;
                    return ERROR_ARG_VALUE;
                }
                argcount+=2;
            }
            else if ( strcasecmp(cpartype,static_cast<const char*>("nrep"))==0 )
            {
                pars->nrep= static_cast<size_t>(std::floor(*mxGetPr(parval)));
                argcount+=2;
            }
            else if ( strcasecmp(cpartype,static_cast<const char*>("seed"))==0 )
            {
                pars->rand_seed = static_cast<int>(std::floor(*mxGetPr(parval)));
                argcount+=2;
            }
            else if ( strcasecmp(cpartype,static_cast<const char*>("verbosity"))==0 )
            {
                pars->verbosity_level = static_cast<int>(std::floor(*mxGetPr(parval)));
                argcount+=2;
            }
            else
            {
                *argposerr = argcount;
                return ERROR_ARG_UNKNOWN;
            }
        }
        else //else return the position of the argument and type of error
        {
            *argposerr = argcount;
            return ERROR_ARG_TYPE;
        }
        mxFree(cpartype); // free the converted argument
    }
    return PACO_NO_ERROR;
}

#include <algorithm>
#include <Eigen/Sparse>
void mexFunction(int nOutputArgs, mxArray *outputArgs[], int nInputArgs, const mxArray * inputArgs[])
{
    PacoParams pars;
    // Set default values for parameters
    pars.qual = QualitySurprise;
    pars.method = MethodAgglomerative;
    pars.nrep = 2; // two are necessary
    pars.verbosity_level=7;
    pars.rand_seed = -1; // default value for the random seed, if -1 then microseconds time is used.

    FILELog::ReportingLevel() = static_cast<TLogLevel>(pars.verbosity_level);

    // Check the arguments of the function
    int error_arg_pos=-1;
    error_type err = parse_args(nOutputArgs, outputArgs, nInputArgs, inputArgs, &pars, &error_arg_pos);

    if (err!=NO_ERROR)
    {
        std::stringstream ss;
        ss << "Error at argument: " << error_arg_pos  << ": " << error_strings[err] ;
        if (err == ERROR_NOT_ENOUGH_ARGS)
            printUsage();
        mexErrMsgTxt(ss.str().c_str());
    }

#ifdef _DEBUG
    printf("[INFO] Method=%d\n[INFO] Consider_comms=%d\n[INFO] CPMgamma=%f\n[INFO] Delta=%f\n[INFO] Max_itr=%zu\n[INFO] Random_order=%d rand_seed=%d\n",pars.method, pars.consider_comms, pars.cpmgamma, pars.delta, pars.max_itr, pars.random_order, pars.rand_seed);
#endif
    // Get number of vertices in the network
    int N = mxGetN(inputArgs[0]); // number of columns
    int M = mxGetM(inputArgs[0]); // number of rows
    double *W = mxGetPr(inputArgs[0]);

    // In this case we are feeding instead of the adjacency matrix, the result of [i j w]=find(A);
    bool feedingSparseMatrix=false;
    if (M>3 && N==3)
    {
        feedingSparseMatrix=true;
    }

    // Create the Graph helper object specifying edge weights too
    GraphC *G=NULL;
    try
    {
        if (feedingSparseMatrix)
        {
            // Map the Mx3 array memory to an Eigen container, to facilitate handling
            Eigen::MatrixXd IJW = Eigen::Map<Eigen::MatrixXd>(W,M,N);
            Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic> B1,B2;
            B1 = (IJW.col(0).array() >= IJW.col(1).array()).cast<int>();
            B2 = (IJW.col(0).array() < IJW.col(1).array()).cast<int>();

            bool isUpperTriangular=false;
            bool isLowerTriangular=false;
            bool isSymmetric=false;
            int sum1 = B1.sum();
            int sum2 = B2.sum();
            // Condizione semplice da verificare facendo [i j w]=find(A), oppure [i j w]=find(triu(A)) oppure [i j w]=find(tril(A))
            if (sum1 == sum2)
                isSymmetric=true;
            if (sum1==0 && sum2==M)
                isUpperTriangular=true;
            if (sum1==M && sum2==0)
                isLowerTriangular=true;

            //printf("sum1=%d sum2=%d Is Symmetric=%d Is isUpperTriangular=%d isLowerTriangular=%d\n",sum1,sum2,isSymmetric,isUpperTriangular,isLowerTriangular);
            
            if (!isSymmetric && !isUpperTriangular && !isLowerTriangular)
            {
                throw std::logic_error("Matrix is not symmetric, nor triangular lower or upper triangular. Check diagonal and non symmetric values.");
            }

            std::vector<double> edges_list;
            std::vector<double> edges_weights;

            for (int l=0; l<M; ++l)
            {
                double row_node = IJW(l,0); //indice riga della find
                double column_node = IJW(l,1); //indice colonna della find
                double w = IJW(l,2);

                if ( isUpperTriangular || isLowerTriangular) // keeps only symmetric and also avoid self-loops (implicitly inserting upper triangular)
                {
                    edges_list.push_back(column_node-1);
                    edges_list.push_back(row_node-1);
                    edges_weights.push_back(w);
                }
                else if (isSymmetric)
                {
                    if (row_node<column_node)
                    {
                        edges_list.push_back(column_node-1);
                        edges_list.push_back(row_node-1);
                        edges_weights.push_back(w);
                    }
                }
                ctrlcCheckPoint(__FILE__, __LINE__); // Interrupt here
            }
            
            G = new GraphC(edges_list.data(),edges_weights.data(),edges_weights.size());
        }
        else
        {
            G  = new GraphC(mxGetPr(inputArgs[0]),N,N);
        }
        // Create an instance of the optimizer
        CommunityStructure c(G);
        c.set_random_seed(pars.rand_seed);
        double finalquality=c.optimize(pars.qual,pars.method,pars.nrep);
        // Prepare output
        outputArgs[0] = mxCreateDoubleMatrix(1,(mwSize)G->number_of_nodes(), mxREAL);
        // Copy the membership vector to outputArgs[0] which has been already preallocated
        igraph_vector_copy_to(c.get_membership(),mxGetPr(outputArgs[0]));
        // Copy the value of partition quality
        outputArgs[1] = mxCreateDoubleScalar(finalquality);
        // Cleanup the memory (follow this order)
        delete G;
    }
    catch (std::exception &e)
    {
        cerr << e.what() << endl;
        mexErrMsgTxt(e.what());
    }

    // Finish the function
    return;
}
