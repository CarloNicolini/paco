#include <iostream>
#include <string.h>
#include <sstream>
#include <igraph.h>
#include <sys/time.h>

// Include MEX for matlab or Octave
#include "mex.h"

#include "Community.h"
#include "Graph.h"
#include "RandomOptimizer.h"
#include "AgglomerativeOptimizer.h"
#include "AnnealOptimizer.h"
#include "SurpriseFunction.h"
#include "AsymptoticSurpriseFunction.h"
#include "SignificanceFunction.h"
#include "igraph_utils.h"
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
    "Non valid input adjacency matrix. Must be symmetric real square matrix.",
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

    if (M!=N  || mxIsComplex(W) || mxIsEmpty(W) || mxIsCell(W) || !mxIsNumeric(W))
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

            mexPrintf("ARGUMENT: %s VALUE=%g\n", cpartype,*mxGetPr(parval));
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

int main(int argc, char *argv[])
{
    return 0;
}
