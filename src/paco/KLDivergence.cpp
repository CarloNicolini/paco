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

#include <cmath>

/**
 * @brief KLs Binary Kullback Leibler divergence between 
 * two binomial distributions with parameters p and q respectively
 * @param q parameter of left binomial distribution
 * @param p parameter of right binomial distribution
 * @return the KL divergence
 */
double KL(double q, double p)
{
    double KL = 0.0;
    if (q > 0.0 && p > 0.0)
        KL += q*log10(q/p);
    if (q < 1.0 && p < 1.0)
        KL += (1.0-q)*log10((1.0-q)/(1.0-p));
    return KL;
}
