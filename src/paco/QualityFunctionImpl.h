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


#ifndef QUALITYFUNCTIONIMPL_H
#define QUALITYFUNCTIONIMPL_H

#include <iostream>

#include <igraph.h>
#include "igraph_utils.h"

#include "PartitionHelper.h"

using std::cout;
using std::cerr;
using std::endl;

class QualityFunctionImpl
{
public:

    virtual ~QualityFunctionImpl() {}

    double &operator()()
    {
        return quality;
    }

    double &operator()(const igraph_t *g, const igraph_vector_t *m, const igraph_vector_t *weights=NULL) const
    {
        eval(g,m,weights);
        return quality;
    }

    double &operator()(const PartitionHelper *par) const
    {
        eval(par);
        return quality;
    }

    virtual QualityFunctionImpl* clone() const = 0;

//protected:
    mutable double quality; // mutable keyword allows to modify quality even in const methods.
    virtual void eval(const igraph_t *g, const igraph_vector_t *m, const igraph_vector_t *weights=NULL) const = 0;
    virtual void eval(const PartitionHelper *par) const = 0;
    PartitionHelper *par;

};

#endif // QUALITYFUNCTIONIMPL_H

