/* child process*/
#pragma comment (lib, "msmpi.lib")
#define MPICH_SKIP_MPICXX
#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>
#include "variable.h"



using namespace std;

int main(int argc, char *argv[])
{
	int rank, i = 0, j = 0, MAXGEN = 0, iter = 0, update = 0, flag = 0;
	int count = 0;
	double nointeg = 0;
	double uninteg = 0.0;
	double inputs[Multipurpose][vardim] = { 0 };  //假定每一类中测试数据的个数,该值应足够大
	double outputs[Multipurpose] = { 0 };      //假定每一类中适应值个数,该值应等于测试数据个数
	double cen[vardim] = { 0 };
	double var[vardim] = { 0 };
	double sum = 0.0;
	double width = 0.0;
	double w = 0.0;
	double out1 = 0.0;
	vector<double> out2;
	double error = 10000;
	double weight = 0.0;
	vector<double> estfit;
	double pop[2 * popsize][vardim] = { 0 };
	double popfit[2 * popsize] = { 0 };
	double person[vardim] = { 0 };
	double fitness = 0.0;
	MPI_Request req;
	MPI_Status status;
	MPI_Comm icomm2;
	MPI_Init(&argc, &argv);
	MPI_Comm_get_parent(&icomm2);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);


	if (rank == 0)
	{
		MPI_Recv(&MAXGEN, 1, MPI_INT, MASTER, 22, icomm2, &status);
	}
	MPI_Bcast(&MAXGEN, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	for (iter = 0; iter <MAXGEN; iter++)
	{
		if (rank == 0)
		{
			MPI_Recv(&update, 1, MPI_INT, MASTER, iter, icomm2, &status);
		}
		MPI_Bcast(&update, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		if (update == 1)
		{
			MPI_Recv(inputs, Multipurpose * vardim, MPI_DOUBLE, MASTER, rank + 1, icomm2, &status);
			MPI_Recv(outputs, Multipurpose, MPI_DOUBLE, MASTER, rank + 2, icomm2, &status);
			count = 0;
			for (i = 0; i < Multipurpose; i++)
			{
				nointeg = 0.0;
				for (j = 0; j < vardim; j++)
				{
					nointeg += inputs[i][j];
				}
				if (nointeg > thres)
				{
					count++;
				}
			}
			if (count == 0)
			{
				count = 1;
			}
			for (i = 0; i < vardim; i++)
			{
				uninteg = 0.0;
				for (j = 0; j < count; j++)
				{
					uninteg += inputs[j][i];
				}
				cen[i] = uninteg / count;
				//cout << "cen[" << i << "]" << cen[i] << endl;
			}
			sum = 0.0;
			for (i = 0; i < count; i++)
			{
				uninteg = 0.0;
				for (j = 0; j < vardim; j++)
				{
					var[j] = pow((inputs[i][j] - cen[j]), 2);
				}
				for (j = 0; j < vardim; j++)
				{
					uninteg += var[j];
				}
				sum += uninteg;
			}
			width = sum / count;
			error = 10000;
			for (w = 0; w - 10 < 0.01; w += 0.1)
			{
				out2.clear();
				for (i = 0; i < count; i++)
				{
					uninteg = 0.0;
					for (j = 0; j < vardim; j++)
					{
						uninteg += pow((inputs[i][j] - cen[j]), 2);
					}
					out1 = exp(-uninteg / (2 * width));
					out2.push_back(w*out1);
				}
				uninteg = 0.0;
				for (i = 0; i < count; i++)
				{
					uninteg += pow((outputs[i] - out2[i]), 2);
				}
				if (uninteg < error)
				{
					weight = w;
					error = uninteg;
					//cout << "weight=" << weight << endl;
				}
			}
			//////////////////////////////////////////////////////////////////////////////
			MPI_Recv(pop, 2 * popsize * vardim, MPI_DOUBLE, MASTER, rank + 3, icomm2, &status);
			count = 0;
			for (i = 0; i < 2 * popsize; i++)
			{
				nointeg = 0.0;
				for (j = 0; j < vardim; j++)
				{
					nointeg += pop[i][j];
				}
				if (nointeg > thres)
				{
					count++;
				}
			}
			estfit.clear();
			for (i = 0; i < count; i++)
			{
				uninteg = 0.0;
				for (j = 0; j < vardim; j++)
				{
					uninteg += pow((pop[i][j] - cen[j]), 2);
				}
				out1 = exp(-uninteg / (2 * width));
				estfit.push_back(weight*out1);
			}
			for (i = 0; i < count; i++)
			{
				popfit[i] = estfit[i];
				//cout << "popfit[" << i << "]" << popfit[i] << endl;
			}
			MPI_Send(popfit, 2 * popsize, MPI_DOUBLE, MASTER, rank + 4, icomm2);
			/////////////////////////////以上是接收属于该类的进化个体///////////////////////////
		}
		else
		{
			MPI_Recv(pop, 2 * popsize * vardim, MPI_DOUBLE, MASTER, rank + 3, icomm2, &status);
			count = 0;
			for (i = 0; i < 2 * popsize; i++)
			{
				nointeg = 0.0;
				for (j = 0; j < vardim; j++)
				{
					nointeg += pop[i][j];
				}
				if (nointeg > thres)
				{
					count++;
				}
			}
			estfit.clear();
			for (i = 0; i < count; i++)
			{
				uninteg = 0.0;
				for (j = 0; j < vardim; j++)
				{
					uninteg += pow((pop[i][j] - cen[j]), 2);
				}
				out1 = exp(-uninteg / (2 * width));
				estfit.push_back(weight*out1);
			}
			for (i = 0; i < count; i++)
			{
				popfit[i] = estfit[i];
				//cout << "popfit[" << i << "]" << popfit[i] << endl;
			}
			MPI_Send(popfit, 2 * popsize, MPI_DOUBLE, MASTER, rank + 4, icomm2);
			/////////////////////////////以上是接收属于该类的进化个体///////////////////////////
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if (rank == 0)
		{
			MPI_Recv(&flag, 1, MPI_INT, MASTER, 33, icomm2, &status);
		}
		MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (flag == 1)
			break;
	}
	// 在主程序结束时，断开通信器
	if (icomm2 != MPI_COMM_NULL) {
		MPI_Barrier(MPI_COMM_WORLD);       // 进行全局同步
		MPI_Barrier(icomm2);               // 父进程同步
		MPI_Comm_disconnect(&icomm2);      // 断开主通信器
		icomm2 = MPI_COMM_NULL;            // 清空句柄
	}
	MPI_Finalize();
	return 0;
}

/*
if (update == 0)
{

}
*/