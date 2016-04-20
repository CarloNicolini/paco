#include <iostream>
#include "Surprise.h"
using namespace std;

int main(int argc, char *argv[])
{
	long int p=atoi(argv[1]);
	long int pi=atoi(argv[2]);
	long int m=atoi(argv[3]);
	long int mi=atoi(argv[4]);
	long int ni=atoi(argv[5]);
		
	double S = computeConditionedSurprise(p,pi,m,mi,ni);
	std::cout << S << std::endl;
	return 0;
}