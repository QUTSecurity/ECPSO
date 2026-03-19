#pragma comment (lib, "msmpi.lib")
#define  MPICH_SKIP_MPICXX 
#include <iostream>
#include <mpi.h>
#include <cmath>
#include <vector>
#include <windows.h>
#include "variable.h"
using namespace std;



void proposedalgorithm(struct individual populationALL[], double sampling[sampsize + selsize*MAXGEN][vardim], double outputing[sampsize + selsize*MAXGEN], struct  individual population[], double sampling2[sampsize + selsize*MAXGEN][vardim], double outputing2[sampsize + selsize*MAXGEN], struct  individual population2[], int update, MPI_Comm icomm2, double a)
{
	int count = 0;
	int i = 0, j = 0, k = 0, l = 0, m = 0, n = 0.0, w = 0.0;
	double tempdouble1 = 0.0;
	double tempdouble2 = 0.0;
	vector<double> squares;
	vector<double> dis_bou_xmin;
	vector<double> dis_ind_xmin;
	vector<double> dis_cen_ind;
	vector<double> dis_samp_xmin;
	vector<double> weight;
	vector<double> fit;
	double inputs[clustersnum][vardim] = { 0 };       
	double outputs[clustersnum] = { 0 };           
	double pop[2 * popsize][vardim] = { 0 };
	double popfit[2 * popsize] = { 0 };
	double pop2[2 * popsize][vardim] = { 0 };
	double popfit2[2 * popsize] = { 0 };
	double tempop[popsize][vardim] = { 0 };
	double tempopfit[popsize] = { 0 };
	double tempop2[popsize][vardim] = { 0 };
	double tempopfit2[popsize] = { 0 };
	MPI_Status status;

	
	if (update == 1)
	{
		MPI_Send(sampling, clustersnum * vardim, MPI_DOUBLE, 0, 0 + 1, icomm2);
		MPI_Send(outputing, clustersnum, MPI_DOUBLE, 0, 0 + 2, icomm2);
		MPI_Send(sampling2, clustersnum * vardim, MPI_DOUBLE, 1, 1 + 1, icomm2);
		MPI_Send(outputing2, clustersnum, MPI_DOUBLE, 1, 1 + 2, icomm2);
	}
	l = 0;


	
	memset(pop, 0, sizeof(pop));
	memset(popfit, 0, sizeof(popfit));
	memset(pop2, 0, sizeof(pop2));
	memset(popfit2, 0, sizeof(popfit2));
	for (j = 0; j < popsize; j++)
	{
		for (k = 0; k < vardim; k++)
		{
			pop[l][k] = population[j].p[k];
			pop2[l][k] = population[j].p[k];
		}
		l++;
	}
	for (j = 0; j < popsize; j++)
	{
		for (k = 0; k < vardim; k++)
		{
			pop[l][k] = population2[j].p[k];
			pop2[l][k] = population2[j].p[k];
		}
		l++;
	}
	MPI_Send(pop, 2 * popsize * vardim, MPI_DOUBLE, 0, 0 + 3, icomm2);
	MPI_Recv(popfit, 2 * popsize, MPI_DOUBLE, 0, 0 + 4, icomm2, &status);
	MPI_Send(pop2, 2 * popsize * vardim, MPI_DOUBLE, 1, 1 + 3, icomm2);
	MPI_Recv(popfit2, 2 * popsize, MPI_DOUBLE, 1, 1 + 4, icomm2, &status);
	
	
	for (j = 0; j < popsize; j++)
	{
		tempdouble1 = 0.0;
		tempdouble2 = 0.0;
		for ( k = 0; k < vardim; k++)
		{
			tempdouble1 += fabs(pop[j][k]);
			tempdouble2 += fabs(pop[j + popsize][k]);
		}

		
		if (tempdouble1 > thres)
		{
			for (k = 0; k < vardim; k++)
			{
				tempop[m][k] = pop[j][k];
			}
			tempopfit[m] = a * popfit[j] + (1 - a) * popfit2[j];
			m++;
		}
		if (tempdouble2 > thres)
		{
			for (k = 0; k < vardim; k++)
			{
				tempop2[n][k] = pop[j + popsize][k];
			}
			tempopfit2[n] = a * popfit[j + popsize] + (1 - a) * popfit2[j + popsize];
			n++;
		}
	}
	w = 0;
	for (i = 0; i < popsize; i++)
	{
		for (j = 0; j < vardim; j++)
		{
			population[i].p[j] = tempop[i][j];
			population2[i].p[j] = tempop2[i][j];
			populationALL[w].p[j] = tempop[i][j];
		}
		population[i].fitness = tempopfit[i];
		population2[i].fitness = tempopfit2[i];
		populationALL[w].fitness = tempopfit[i];
		w++;
	}
	for (i = 0; i < popsize; i++)
	{
		for (j = 0; j < vardim; j++)
		{
			populationALL[w].p[j] = tempop2[i][j];
		}
		populationALL[w].fitness = tempopfit2[i];
		w++;
	} 
}
