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


#ifndef _ANNEALOPTIMIZER_H
#define _ANNEALOPTIMIZER_H

#include "QualityOptimizer.h"

class AnnealParameters
{
public:
    AnnealParameters(size_t nIter=1E3, size_t maxHits=1E2, double mintemp=1E-1, double tol=1E-5, double temp=1E6, double tempscale=0.99, double minfval=std::numeric_limits<double>::max()) :
        nIterations(nIter),
        nHits(maxHits),
        min_temp(mintemp),
        tolerance(tol),
        temperature(temp),
        temp_scale(tempscale),
        minfval(minfval)
    {}
    size_t nIterations;
    size_t nHits;
    double min_temp;
    double tolerance;
    double temperature;
    double temp_scale;
    double minfval;
};

class AnnealOptimizer : public QualityOptimizer
{
public:
    AnnealOptimizer() {}
    AnnealOptimizer(const igraph_t *g, const QualityFunction &fun, const igraph_vector_t *memb, const igraph_vector_t *weights=NULL);
    virtual ~AnnealOptimizer();
    void set_parameters(const AnnealParameters &_par)
    {
        this->param = _par;
    }
    double optimize(const igraph_t *g, const QualityFunction &fun, const igraph_vector_t *memb, const igraph_vector_t *weights=NULL);

protected:
    double diff_move(const igraph_t *g, const QualityFunction &fun, const igraph_vector_t *memb, int vert, size_t dest_comm, const igraph_vector_t *weights);

    AnnealParameters param;
};
#endif // _ANNEALOPTIMIZER_H
