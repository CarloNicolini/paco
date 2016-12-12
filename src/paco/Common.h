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


#ifndef GRAPH_COMMON
#define GRAPH_COMMON

#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <unordered_map>
#include <unordered_set>


// http://stackoverflow.com/questions/3219393/stdlib-and-colored-output-in-c
#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

using std::cout;
using std::endl;
using std::cerr;

using std::vector;
using std::map;
using std::unordered_map;
using std::unordered_set;
using std::ostream;
using std::pair;
using std::set;
using std::string;
using std::ifstream;
using std::ofstream;

/**
 * @brief num_pairs
 * @param x
 * @return
 */
inline int num_pairs(int x)
{
    if (x<=1)
        return 0;
    else
        return x*(x-1)/2;
}


/**
 * @brief mapvalue_sum
 * @param m
 * @return
 */
template <class T>
T mapvalue_sum(const map<size_t,T> &m)
{
    T sum = 0;
    for ( typename map<size_t,T>::const_iterator i = m.begin(); i!=m.end(); ++i)
        sum+= i->second;
    return sum;
}

template <class T>
T mapvalue_sum(const unordered_map<size_t,T> &m)
{
    T sum = 0;
    for ( typename unordered_map<size_t,T>::const_iterator i = m.begin(); i!=m.end(); ++i)
        sum+= i->second;
    return sum;
}

/**
 * @brief operator <<
 * @param os
 * @param m
 * @return
 */
template <class T>
inline ostream &operator<<(ostream &os, const map<size_t,T> &m)
{
    for ( typename map<size_t,T>::const_iterator i=m.begin() ; i!= m.end() ; ++i)
    {
        os << i->first << "-> " << i->second << endl;
    }
    return os;
}


/**
 * @brief operator <<
 * @param os
 * @param x
 * @return
 */
template <class T>
inline ostream& operator<<(ostream& os, vector<T> x)
{
    os << "[ ";
    for ( typename vector<T>::const_iterator it = x.begin(); it!=x.end(); ++it)
        os << *it << ", ";
    os << "]" << endl;
    return os;
}

/**
 * @brief sort_second
 * @param l
 * @param r
 * @return
 */
template <class T>
inline bool sort_second( const pair<int,T>& l,  const pair<int,T>& r)
{
    return l.second > r.second;
}

#endif // COMMON

