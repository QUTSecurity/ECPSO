#include "variable.h"


float Sim(const int BL[numprocs][CN], const int MB[numprocs][CN])
{
	float subsim[numprocs] = { 0.0 },Fit=0.0;
	int BLlength[numprocs] = { 0 }, MBlength[numprocs] = { 0 }, Max[numprocs] = { 0 };
	int i,j;
	for (i = 0; i < numprocs;i++)
	for (j = 0; j < CN; j++)
	{
		if (BL[i][j] != 0)
			BLlength[i] += 1;
	}
	for (i = 0; i < numprocs; i++)
	for (j = 0; j < CN; j++)
	{ 
		if (MB[i][j] != 0)
			MBlength[i] += 1;
	}
	for (i = 0; i < numprocs; i++)
	{
		if (BLlength[i] < MBlength[i])
			Max[i] = MBlength[i];
		else
			Max[i] = BLlength[i];
	}
	for (i = 0; i < numprocs; i++)
	{
		for (j = 0; j < Max[i]; j++)
		{
			if (BL[i][j] != MB[i][j])
				break;
		}
		subsim[i] = j / (float(Max[i])); 
	}
	for (i = 0; i < numprocs; i++)
	{
		Fit += subsim[i];
	}
	return Fit / numprocs;
}
