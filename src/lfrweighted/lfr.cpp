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

#include "set_parameters.h" // the LFR parameters
#include "../paco/igraph_utils.h" // to handle conversion from EigenMatrix to igraph object and then to mxArray

#ifdef __linux__
    #include <mex.h>
#endif

#ifdef __APPLE__
#include "mex.h"
#endif

using namespace std;

void printUsage()
{
    mexPrintf("LFRW: Lancichinetti-Fortunato-Radicchi network generator\n");

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

/**
 * @brief parse_args
 * @param nOutputArgs
 * @param outputArgs
 * @param nInputArgs
 * @param inputArgs
 * @param pars
 * @param argposerr
 * @return
 */
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

/*
    command_flags.push_back("-N");			//0
    command_flags.push_back("-k");			//1
    command_flags.push_back("-maxk");		//2
    command_flags.push_back("-mut");		//3
    command_flags.push_back("-t1");			//4
    command_flags.push_back("-t2");			//5
    command_flags.push_back("-minc");		//6
    command_flags.push_back("-maxc");		//7
    command_flags.push_back("-on");			//8
    command_flags.push_back("-om");			//9
    command_flags.push_back("-beta");		//10
    command_flags.push_back("-muw");		//11
    command_flags.push_back("-C");			//12
*/

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
            if ( strcasecmp(cpartype,"N")==0 )
            {
                /*
                pars->method = static_cast<OptimizerType>(*mxGetPr(parval));
                if (pars->method<0 || pars->method>3)
                {
                    *argposerr = argcount+1;
                    return ERROR_ARG_VALUE;
                }
                */
                argcount+=2;
            }
            else if ( strcasecmp(cpartype,"Quality")==0 )
            {
                /*
                pars->qual = static_cast<QualityType>(*mxGetPr(parval));
                if (pars->qual<0 || pars->qual>3)
                {
                    *argposerr = argcount+1;
                    return ERROR_ARG_VALUE;
                }
                */
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

/**
 * @brief mexFunction
 * @param nOutputArgs
 * @param outputArgs
 * @param nInputArgs
 * @param inputArgs
 */
void mexFunction(int nOutputArgs, mxArray *outputArgs[], int nInputArgs, const mxArray * inputArgs[])
{
    FILELog::ReportingLevel() = static_cast<TLogLevel>(7);
    // Set standard parameters argument

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

    int N=1000;
    try
    {
       // Prepare output
        outputArgs[0] = mxCreateDoubleMatrix((mwSize)N,(mwSize)N, mxREAL);
        // Copy the membership vector to outputArgs[0] which has been already preallocated
        igraph_vector_copy_to(c.get_membership(),mxGetPr(outputArgs[0]));

        // Copy the value of partition quality
        outputArgs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
        double *q = mxGetPr(outputArgs[1]);
        // Copy final quality value
        q[0] = finalquality;
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
