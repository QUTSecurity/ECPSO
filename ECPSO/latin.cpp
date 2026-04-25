#include <iostream>
#include <vector>
#include <random>
#include "variable.h"

using namespace std;
using std::default_random_engine;
using std::uniform_int_distribution;


void latin_hypercube(double  samples[sampsize][vardim])
{ 
	int i, j;
	double d = 1.0 / sampsize;
	double nosie = 0;
	vector< vector <double> > sampling(vardim, vector<double>(sampsize));
	double temp[sampsize];
	default_random_engine e;
	for (j = 0; j < sampsize; j++)
	{
		
		uniform_real_distribution<double> u(j * d + nosie, (j + 1) * d - nosie); 
		temp[j] = u(e);
	}
	if (xmax < xmin)
	{
		cout << "The range is wrong" << endl;
	}
	for (i = 0; i < vardim; i++)
	{
		for (j = 0; j < sampsize; j++)
		{
			uniform_real_distribution<double> u((temp[j] * (xmax - xmin)), ((temp[j] + d)*(xmax - xmin)));
			sampling[i][j] = u(e);
		}
	}
	for (i = 0; i < vardim; i++)
	{
		random_shuffle(sampling[i].begin(), sampling[i].end());
	}
	for (i = 0; i < sampsize; i++)
	{
		for (j = 0; j < vardim; j++)
		{
			samples[i][j] = sampling[j][i];
		}
	}

}
