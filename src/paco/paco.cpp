/* This file is part of PACO-PArtitioning Clustering Optimization a program
* to find network partitions using modular solvers and quality functions.
*
*  Copyright (C) 2015 Carlo Nicolini <carlo.nicolini@iit.it>
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
*/

#include <iostream>
#include <string.h>
#include <sstream>
#include <igraph.h>
#include <sys/time.h>

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
	#ifdef MATLAB_SUPPORT
	#include "/usr/local/MATLAB/R2015a/extern/include/mex.h"
	#elif OCTAVE_SUPPORT
	#include "/usr/include/octave-4.0.0/octave/mex.h"
	#endif
#endif

#ifdef __APPLE__
#include "mex.h"
#endif

using namespace std;

enum PacoMethod
{
    MethodAgglomerative = 0,
    MethodRandom = 1,
    MethodAnneal = 2,
};

enum PacoQuality
{
    QualitySurprise = 0,
    QualitySignificance = 1,
    QualityAsymptoticSurprise = 2
};

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
    mexPrintf("	quality is one of the following integers: {0,1,2}:\n");
    mexPrintf("		0: Surprise (discrete)\n");
    mexPrintf("		1: Significance\n");
    mexPrintf("		2: AsymptoticSurprise\n");
    mexPrintf("[m, qual] = paco(W,'seed',val);\n");
    mexPrintf("		val: to provide a specific random seed to the algorithm, in order to have reproducible results.\n");
    mexPrintf("\n\n");
    mexPrintf("Example:\n");
    mexPrintf(">> A=rand(100,100); A=(A+A')/2; A=A.*(A>0.5);\n");
    mexPrintf(">> [memb, qual] = paco(A,'method',2);\n");
}

enum error_type
{
    NO_ERROR = 0,
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
    "Non valid input adjacency matrix. Must be symmetric real dense-type square matrix.",
    "Expected some argument value but empty found.",
    "Unkwown argument."
};

struct PacoParams
{
    PacoMethod method;
    PacoQuality qual;
    size_t nrep;      // Maximum number of consecutive repetitions to perform.
    int rand_seed; // random seed for the louvain algorithm
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

    bool v1 = M!=N;
    bool v2 = mxIsComplex(W);
    bool v3 = mxIsEmpty(W);
    bool v4 = mxIsCell(W);
    bool v5 = !mxIsNumeric(W);
    bool v6 = !mxIsSparse(W);

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
                pars->method = static_cast<PacoMethod>(*mxGetPr(parval));
                if (pars->method<0 || pars->method>2)
                {
                    *argposerr = argcount+1;
                    return ERROR_ARG_VALUE;
                }
                argcount+=2;
            }
            else if ( strcasecmp(cpartype,"Quality")==0 )
            {
                pars->qual = static_cast<PacoQuality>(*mxGetPr(parval));
                if (pars->qual<0 || pars->qual>2)
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
    return NO_ERROR;
}

void mexFunction(int nOutputArgs, mxArray *outputArgs[], int nInputArgs, const mxArray * inputArgs[])
{
    FILELog::ReportingLevel() = static_cast<TLogLevel>(7);
    PacoParams pars;
    // Set default values for parameters
    pars.qual = QualitySurprise;
    pars.method = MethodAgglomerative;
    pars.nrep = 1;
    pars.rand_seed = -1; // default value for the random seed, if -1 than microseconds time is used.

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
    int N = mxGetN(inputArgs[0]);

    // Create the Graph helper object specifying edge weights too
    GraphC *G;
    try
    {
        G  = new GraphC(mxGetPr(inputArgs[0]),N,N);
    }
    catch (std::exception &e)
    {
        //mexErrMsgTxt("Input network has diagonal entries. Set them to zero.");
        std::string error_string = e.what();
        delete G;
        mexErrMsgTxt(error_string.c_str());
    }

    // Create an instance of the optimizer
    CommunityStructure c(G);
    c.set_random_seed(pars.rand_seed);
    c.sort_edges();

    QualityFunction *fun;
    QualityOptimizer *opt;

    // Select the optimization method
    switch (pars.method)
    {
    case MethodAgglomerative:
    {
        opt = dynamic_cast<AgglomerativeOptimizer*>( new AgglomerativeOptimizer);
        dynamic_cast<AgglomerativeOptimizer*>(opt)->set_edges_order(c.get_sorted_edges_indices());
        break;
    }
    case MethodRandom:
    {
        opt = dynamic_cast<RandomOptimizer*>(new RandomOptimizer);
        break;
    }
    case MethodAnneal:
    {
        opt = dynamic_cast<AnnealOptimizer*>(new AnnealOptimizer);
        break;
    }
    default:
    {
        mexErrMsgTxt("Non supported optimization method");
        break;
    }
    }

    // Select the quality function
    switch (pars.qual)
    {
    case QualitySurprise:
    {
        fun = dynamic_cast<SurpriseFunction*>(new SurpriseFunction);
        break;
    }
    case QualitySignificance:
    {
        fun = dynamic_cast<SignificanceFunction*>(new SignificanceFunction);
        break;
    }
    case QualityAsymptoticSurprise:
    {
        fun = dynamic_cast<AsymptoticSurpriseFunction*>(new AsymptoticSurpriseFunction);
        break;
    }
    default:
    {
        mexErrMsgTxt("Non supported quality function");
        break;
    }
    }

    // Finally optimize the partition
    bool is_weighted = G->is_weighted();
    if (is_weighted && pars.qual == QualitySurprise)
        mexErrMsgTxt("Can't optimize discrete surprise on weighted graph. Use AsymptoticSurprise instead.");

    double finalqual = 0;
    try
    {
        for (int i=0; i<pars.nrep; ++i)
        {
            if (!is_weighted)
                finalqual = opt->optimize(G->get_igraph(),*fun,c.get_membership());
            else
                finalqual = opt->optimize(G->get_igraph(),*fun,c.get_membership(),G->get_edge_weights());
        }
    }
    catch(std::exception &e)
    {
        mexErrMsgTxt(e.what());
    }

    // Prepare output
    outputArgs[0] = mxCreateDoubleMatrix(1,(mwSize)N, mxREAL);
    double *memb = mxGetPr(outputArgs[0]);

    // Copy the membership vector to outputArgs[0] which has been already preallocated
    igraph_vector_copy_to(c.get_membership(),memb);
    //for (size_t i = 0; i<N ; ++i)
    //  memb[i] = static_cast<double>(c.get_membership(i));

    // Copy the value of partition quality
    outputArgs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
    double *q = mxGetPr(outputArgs[1]);

    // Copy final quality value
    q[0] = finalqual;


    // Cleanup the memory (follow this order)
    delete G;
    delete opt;
    delete fun;
    // Finish the function
    return;
}
