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

#include <Eigen/Core>

using namespace std;


int main()
{
  Eigen::MatrixXd W;
  W.setRandom(100,100);
  for (int i=0; i<100; ++i)
    W(i,i)=0;
  //W = W+W.transpose();
  //W.cwiseSum(2.0);
  GraphC *G;
  G = new GraphC(W);
  delete G;
  return 0;
}
