#pragma comment (lib, "msmpi.lib")
#define  MPICH_SKIP_MPICXX 
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "variable.h"
using namespace std;
 
int my_rank, group_size, sub_group;
int n_all;
double pointx[30], pointy[30];
double tempx[2], tempy[2], xtemp[4], ytemp[4];
MPI_Status status;
int location[30];
int tt[CN] = { 0 };
int BL[numpath][CN] = { 0 };
int k = 0;


void getymin()
{
	int i, index = 0;
	double tempvalue = pointy[0];
	tt[k++] = 111;

	for (i = 1; i<n_all; i++)
	{
		tt[k++] = 112;
		tt[k++] = 113;
		if (pointy[i]<tempvalue)
		{
			tempvalue = pointy[i];
			index = i;
			tt[k++] = 114;
		}
		else
		{
			
			tt[k++] = 115;
			if (pointy[i] == tempvalue)
			{
				tt[k++] = 116;
				if (pointx[index] > pointx[i])
				{
					index = i;
					tt[k++] = 117;
				}
			}
		}
	}
	tempx[0] = pointx[index];
	tempy[0] = pointy[index];

	
	tempx[1] = tempx[0];
	tempy[1] = tempy[0];
	tt[k++] = 118;

	
	for (i = 0; i<n_all; i++)
	{
		tt[k++] = 119;
		tt[k++] = 120;
		if ((pointy[i] == tempy[1]) && (pointx[i]>tempx[1]))
		{
			
			tempx[1] = pointx[i];
			tempy[1] = pointy[i];
			tt[k++] = 121;
		}
	}
	return;
}



void getymax()
{
	int i, index = 0;
	double temp = pointy[0];
	tt[k++] = 122;

	for (i = 1; i<n_all; i++)
	{
		tt[k++] = 123;
		tt[k++] = 124;
		if (pointy[i]>temp)
		{
			temp = pointy[i];
			index = i;
			tt[k++] = 125;
		}
		else
		{
			
			tt[k++] = 126;
			if (temp == pointy[i])
			{
				tt[k++] = 127;
				if (pointx[i] > pointx[index])
				{
					index = i;
					tt[k++] = 128;
				}
			}

		}
	}
	
	tempx[0] = pointx[index];
	tempy[0] = pointy[index];

	
	tempx[1] = tempx[0];
	tempy[1] = tempy[0];
	tt[k++] = 129;

	
	for (i = 0; i<n_all; i++)
	{
		tt[k++] = 130;
		tt[k++] = 131;
		if ((pointy[i] == tempy[1]) && (pointx[i]<tempx[1]))
		{
			
			tempx[1] = pointx[i];
			tempy[1] = pointy[i];
			tt[k++] = 132;
		}
	}

	return;
}



void getxmin()
{
	int i, index = 0;
	double temp = pointx[0];
	tt[k++] = 133;

	for (i = 1; i<n_all; i++)
	{
		tt[k++] = 134;
		tt[k++] = 135;
		if (pointx[i]<temp)
		{
			temp = pointx[i];
			index = i;
			tt[k++] = 136;
		}
		else
		{
			tt[k++] = 137;
			
			if (pointx[i] == temp)
			{
				tt[k++] = 138;
				if (pointy[index] < pointy[i])
				{
					index = i;
					tt[k++] = 139;
				}
			}
		}
	}
	
	tempx[0] = pointx[index];
	tempy[0] = pointy[index];

	
	tempx[1] = tempx[0];
	tempy[1] = tempy[0];
	tt[k++] = 140;

	
	for (i = 0; i<n_all; i++)
	{
		tt[k++] = 141;
		tt[k++] = 142;
		if ((pointx[i] == tempx[1]) && (pointy[i]<tempy[1]))
		{
			
			tempx[1] = pointx[i];
			tempy[1] = pointy[i];
			tt[k++] = 143;
		}
	}
	return;
}



void getxmax()
{
	int i, index = 0;
	double temp = pointx[0];
	tt[k++] = 144;

	for (i = 1; i<n_all; i++)
	{
		tt[k++] = 145;
		tt[k++] = 146;
		if (pointx[i]>temp)
		{
			temp = pointx[i];
			index = i;
			tt[k++] = 147;
		}
		else
		{
			tt[k++] = 148;
			
			if (temp == pointx[i])
			{
				tt[k++] = 149;
				if (pointy[index] > pointy[i])
				{
					index = i;
					tt[k++] = 150;
				}
			}
		}
	}
	
	tempx[0] = pointx[index];
	tempy[0] = pointy[index];

	
	tempx[1] = tempx[0];
	tempy[1] = tempy[0];
	tt[k++] = 151;

	
	for (i = 0; i<n_all; i++)
	{
		tt[k++] = 152;
		tt[k++] = 153;
		if ((pointx[i] == tempx[1]) && (pointy[i]>tempy[1]))
		{
			
			tempx[1] = pointx[i];
			tempy[1] = pointy[i];
			tt[k++] = 154;
		}
	}
	return;
}





void getincludedvertex(int tag)
{
	double linea, lineb, dist;
	int i, count = 0;
	tt[k++] = 155;

	
	tt[k++] = 156;
	if ((tempx[1] == tempx[0]) && (tempy[1] == tempy[0]))
	{
		n_all = 1;
		pointx[0] = tempx[1];
		pointy[0] = tempy[1];
		tt[k++] = 157;
	}
	else
	{
		
		linea = (tempy[1] - tempy[0]) / (tempx[1] - tempx[0]);
		lineb = tempy[1] - linea*tempx[1];
		tt[k++] = 158;
		
		for (i = 0; i<n_all; i++)
		{
			tt[k++] = 159;
			dist = linea*pointx[i] + lineb - pointy[i];
			tt[k++] = 160;
			if (tag*dist<0)
			{
				pointx[count] = pointx[i];
				pointy[count++] = pointy[i];
				tt[k++] = 161;
			}
		}
		
		for (i = 0; i<count; i++)
		{
			tt[k++] = 162;
			pointx[count - i] = pointx[count - 1 - i];
			pointy[count - i] = pointy[count - 1 - i];
			tt[k++] = 163;
		}
		
		pointx[0] = tempx[0];
		pointy[0] = tempy[0];
		pointx[count + 1] = tempx[1];
		pointy[count + 1] = tempy[1];
		count += 2;
		
		n_all = count;
		tt[k++] = 164;
	}
}



void nextindex(int i, int x, int y)
{
	double x1, y1, temp, valuemax = 0;
	int j;
	tt[k++] = 165;

	x1 = pointx[i];
	y1 = pointy[i];

	
	for (j = 0; j<n_all; j++)
	{
		tt[k++] = 166;
		tt[k++] = 167;
		if (((pointx[j] - x1)*x + (pointy[j] - y1)*y)>0)
		{
			
			temp = ((pointx[j] - x1)*x + (pointy[j] - y1)*y) / (sqrt((pointx[j] - x1)*(pointx[j] - x1) + (pointy[j] - y1)*(pointy[j] - y1)));
			tt[k++] = 168;
			tt[k++] = 169;
			if (temp>valuemax)
			{
				location[i] = j;
				valuemax = temp;
				tt[k++] = 170;
			}
		}
	}
	tt[k++] = 171;
	if (valuemax == 0)
	{
		location[i] = -1;
		tt[k++] = 172;
	}
}



void output(int rank)
{
	int j = 0, index, temp;
	double xtem, ytem;
	int flag = 0;
	tt[k++] = 173;

	
	tt[k++] = 174;
	if (rank == 0)
	{
		flag = 1;
		tt[k++] = 175;
	}
	tt[k++] = 176;
	if (flag == 0)
	{
		tt[k++] = 177;
		if (!((pointx[n_all - 1] == tempx[0]) && (pointy[n_all - 1] == tempy[0])))
		{
			flag = 1;
			tt[k++] = 178;
		}
	}

	
	index = location[j];
	tt[k++] = 179;
	while (location[index] != -1)
	{
		tt[k++] = 180;
		temp = location[index];
		location[index] = j;
		j = index;
		index = temp;
		tt[k++] = 181;
	}
	location[index] = j;
	location[0] = -1;
	tt[k++] = 182;

	
	j = n_all - 1;
	while (location[j] != -1)
	{
		tt[k++] = 183;
		tt[k++] = 184;
		if (!((j == n_all - 1) && (flag == 0)))
		{
			
			xtem = pointx[j];
			ytem = pointy[j];
			tt[k++] = 185;
		}
		j = location[j];
		tt[k++] = 186;
	}

	
	tt[k++] = 187;
	if (!((j == n_all - 1) && (flag == 0)))
	{
		
		xtem = pointx[j];
		ytem = pointy[j];
		tt[k++] = 188;
	}

	tt[k++] = 189;
	if ((rank == 3 * sub_group) && !((xtem == tempx[1]) && (ytem == tempy[1])))
		tt[k++] = 190;
	
	return;
}


int main(int argc, char *argv[])
{
	int i, j, flag = 0;
	int n;             
	int looptimes = 0, iter = 0;
	int points[vardim] = { 0 };
	MPI_Comm icomm1;
	MPI_Init(&argc, &argv);
	MPI_Comm_get_parent(&icomm1);
	MPI_Comm_size(MPI_COMM_WORLD, &group_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	
	sub_group = group_size / 4;
	
	if (my_rank == 0)
	{
		MPI_Recv(&looptimes, 1, MPI_INT, MASTER, 11, icomm1, &status);
	}
	MPI_Bcast(&looptimes, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	
	for (iter = 0; iter < looptimes; iter++)
	{
		memset(tt, 0, sizeof(tt));
		memset(BL, 0, sizeof(BL));
		k = 0;
		if (my_rank == 0)
		{
			n_all = vardim;
			MPI_Recv(points, vardim, MPI_INT, MASTER, iter, icomm1, &status);
			for (i = 0; i<(vardim / 2); i++)
			{
				pointx[i] = points[i * 2];
				pointy[i] = points[(i * 2) + 1];
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		
		tt[k++] = 1;
		MPI_Bcast(&n_all, 1, MPI_INT, 0, MPI_COMM_WORLD);
		tt[k++] = 2;
		
		
		tt[k++] = 3;
		if (my_rank == 0)
		{
			for (i = 1; i<4; i++)
			{
				tt[k++] = 4;
				MPI_Send(pointx, n_all, MPI_DOUBLE, i*sub_group, i, MPI_COMM_WORLD);
				tt[k++] = 5;
				MPI_Send(pointy, n_all, MPI_DOUBLE, i*sub_group, i, MPI_COMM_WORLD);
				tt[k++] = 6;
			}
		}
		tt[k++] = 7;
		if ((my_rank>0) && (my_rank<4 * sub_group) && (my_rank%sub_group == 0))
		{
			MPI_Recv(pointx, n_all, MPI_DOUBLE, 0, my_rank / sub_group, MPI_COMM_WORLD, &status);
			tt[k++] = 8;
			MPI_Recv(pointy, n_all, MPI_DOUBLE, 0, my_rank / sub_group, MPI_COMM_WORLD, &status);
			tt[k++] = 9;
		}

		MPI_Barrier(MPI_COMM_WORLD);
		tt[k++] = 10;

		
		
		tt[k++] = 11;
		if (my_rank == 0)
		{
			getymin();
			tt[k++] = 12;
		}
		
		tt[k++] = 13;
		if (my_rank == sub_group)
		{
			getxmax();
			tt[k++] = 14;
		}
		
		tt[k++] = 15;
		if (my_rank == 2 * sub_group)
		{
			getymax();
			tt[k++] = 16;
		}
		
		tt[k++] = 17;
		if (my_rank == 3 * sub_group)
		{
			getxmin();
			tt[k++] = 18;
		}

		MPI_Barrier(MPI_COMM_WORLD);
		tt[k++] = 19;

		
		tt[k++] = 20;
		if ((my_rank>0) && (my_rank<4 * sub_group) && (my_rank%sub_group == 0))
		{
			
			MPI_Send(&tempx[1], 1, MPI_DOUBLE, 0, my_rank, MPI_COMM_WORLD);
			tt[k++] = 21;
			MPI_Send(&tempy[1], 1, MPI_DOUBLE, 0, my_rank, MPI_COMM_WORLD);
			tt[k++] = 22;
		}
		tt[k++] = 23;
		if (my_rank == 0)
		{
			tt[k++] = 24;
			for (i = 1; i<4; i++)
			{
				tt[k++] = 25;
				MPI_Recv(&xtemp[i], 1, MPI_DOUBLE, i*sub_group, i*sub_group, MPI_COMM_WORLD, &status);
				tt[k++] = 26;
				MPI_Recv(&ytemp[i], 1, MPI_DOUBLE, i*sub_group, i*sub_group, MPI_COMM_WORLD, &status);
				tt[k++] = 27;
			}
			xtemp[0] = tempx[1];
			ytemp[0] = tempy[1];
			tempx[1] = xtemp[3];
			tempy[1] = ytemp[3];
			tt[k++] = 28;
			for (i = 1; i<4; i++)
			{
				
				tt[k++] = 29;
				MPI_Send(&xtemp[i - 1], 1, MPI_DOUBLE, i*sub_group, i*sub_group, MPI_COMM_WORLD);
				tt[k++] = 30;
				MPI_Send(&ytemp[i - 1], 1, MPI_DOUBLE, i*sub_group, i*sub_group, MPI_COMM_WORLD);
				tt[k++] = 31;
			}
		}
		else
		{
			tt[k++] = 32;
			if ((my_rank<4 * sub_group) && (my_rank%sub_group == 0))
			{
				
				MPI_Recv(&tempx[1], 1, MPI_DOUBLE, 0, my_rank, MPI_COMM_WORLD, &status);
				tt[k++] = 33;
				MPI_Recv(&tempy[1], 1, MPI_DOUBLE, 0, my_rank, MPI_COMM_WORLD, &status);
				tt[k++] = 34;
			}
		}

		MPI_Barrier(MPI_COMM_WORLD);
		tt[k++] = 35;

		
		

		
		
		tt[k++] = 36;
		if (my_rank == 0)
		{
			getincludedvertex(-1);
			tt[k++] = 37;
		}
		tt[k++] = 38;
		if (my_rank == sub_group)
		{
			getincludedvertex(-1);
			tt[k++] = 39;
		}
		tt[k++] = 40;
		if (my_rank == 2 * sub_group)
		{
			getincludedvertex(1);
			tt[k++] = 41;
		}
		tt[k++] = 42;
		if (my_rank == 3 * sub_group)
		{
			getincludedvertex(1);
			tt[k++] = 43;
		}

		MPI_Barrier(MPI_COMM_WORLD);
		tt[k++] = 44;

		tt[k++] = 45;
		if ((my_rank<4 * sub_group) && (my_rank%sub_group == 0))
		{
			for (i = 1; i<sub_group; i++)
			{
				tt[k++] = 46;
				
				MPI_Send(&n_all, 1, MPI_INT, my_rank + i, i, MPI_COMM_WORLD);
				tt[k++] = 47;
			}
		}
		else
		{
			tt[k++] = 48;
			if (my_rank < 4 * sub_group)
			{
				
				MPI_Recv(&n_all, 1, MPI_INT, my_rank / sub_group*sub_group, my_rank%sub_group, MPI_COMM_WORLD, &status);
				tt[k++] = 49;
			}
		}

		MPI_Barrier(MPI_COMM_WORLD);
		tt[k++] = 50;

		tt[k++] = 51;
		if ((my_rank<4 * sub_group) && (my_rank%sub_group == 0))
		{
			for (i = 1; i<sub_group; i++)
			{
				tt[k++] = 52;
				
				MPI_Send(pointx, n_all, MPI_DOUBLE, my_rank + i, i, MPI_COMM_WORLD);
				tt[k++] = 53;
				MPI_Send(pointy, n_all, MPI_DOUBLE, my_rank + i, i, MPI_COMM_WORLD);
				tt[k++] = 54;
			}
		}
		else
		{
			tt[k++] = 55;
			if (my_rank<4 * sub_group)
			{
				MPI_Recv(pointx, n_all, MPI_DOUBLE, my_rank / sub_group*sub_group, my_rank%sub_group, MPI_COMM_WORLD, &status);
				tt[k++] = 56;
				MPI_Recv(pointy, n_all, MPI_DOUBLE, my_rank / sub_group*sub_group, my_rank%sub_group, MPI_COMM_WORLD, &status);
				tt[k++] = 57;
			}
		}

		MPI_Barrier(MPI_COMM_WORLD);
		tt[k++] = 58;

		
		tt[k++] = 59;
		if ((n_all>1) && (my_rank<4 * sub_group))
		{
			for (i = my_rank%sub_group; i<n_all - 1; i += sub_group)
			{
				tt[k++] = 60;
				tt[k++] = 61;
				if (my_rank / sub_group == 0)
				{
					nextindex(i, -1, 0);
					tt[k++] = 62;
				}
					
				tt[k++] = 63;
				if (my_rank / sub_group == 1)
				{
					nextindex(i, 0, -1);
					tt[k++] = 64;
				}
					
				tt[k++] = 65;
				if (my_rank / sub_group == 2)
				{
					nextindex(i, 1, 0);
					tt[k++] = 66;
				}
					
				tt[k++] = 67;
				if (my_rank / sub_group == 3)
				{
					nextindex(i, 0, 1);
					tt[k++] = 68;
				}
					
			}
		}

		MPI_Barrier(MPI_COMM_WORLD);
		tt[k++] = 69;

		
		tt[k++] = 70;
		if ((my_rank<4 * sub_group) && (my_rank%sub_group != 0))
		{
			for (i = my_rank%sub_group; i < n_all - 1; i += sub_group)
			{
				tt[k++] = 71;
				MPI_Send(&location[i], 1, MPI_INT, my_rank / sub_group*sub_group, my_rank, MPI_COMM_WORLD);
				tt[k++] = 72;
			}
				
		}
		else
		{
			tt[k++] = 73;
			if (my_rank < 4 * sub_group)
			{
				for (i = 1; i < n_all - 1; i++)
				{
					tt[k++] = 74;
					tt[k++] = 75;
					if (i%sub_group != 0)
					{
						MPI_Recv(&location[i], 1, MPI_INT, my_rank + i%sub_group, my_rank + i%sub_group, MPI_COMM_WORLD, &status);
						tt[k++] = 76;
					}
				}
			}
		}

		MPI_Barrier(MPI_COMM_WORLD);
		tt[k++] = 77;
		location[n_all - 1] = -1;
		tt[k++] = 78;

		MPI_Barrier(MPI_COMM_WORLD);
		tt[k++] = 79;

		
		tt[k++] = 80;
		if (my_rank == 0)
		{
			
			output(my_rank);
			tt[k++] = 81;
			
			MPI_Send(&pointx[0], 1, MPI_DOUBLE, sub_group, sub_group, MPI_COMM_WORLD);
			tt[k++] = 82;
			MPI_Send(&pointy[0], 1, MPI_DOUBLE, sub_group, sub_group, MPI_COMM_WORLD);
			tt[k++] = 83;
			
			MPI_Send(&pointx[n_all - 1], 1, MPI_DOUBLE, 3 * sub_group, 3 * sub_group, MPI_COMM_WORLD);
			tt[k++] = 84;
			MPI_Send(&pointy[n_all - 1], 1, MPI_DOUBLE, 3 * sub_group, 3 * sub_group, MPI_COMM_WORLD);
			tt[k++] = 85;
		}
		tt[k++] = 86;
		if (my_rank == sub_group)
		{
			
			MPI_Recv(&tempx[0], 1, MPI_DOUBLE, 0, sub_group, MPI_COMM_WORLD, &status);
			tt[k++] = 87;
			MPI_Recv(&tempy[0], 1, MPI_DOUBLE, 0, sub_group, MPI_COMM_WORLD, &status);
			tt[k++] = 88;
			
			MPI_Send(&pointx[0], 1, MPI_DOUBLE, 2 * sub_group, 2 * sub_group, MPI_COMM_WORLD);
			tt[k++] = 89;
			MPI_Send(&pointy[0], 1, MPI_DOUBLE, 2 * sub_group, 2 * sub_group, MPI_COMM_WORLD);
			tt[k++] = 90;
		}
		tt[k++] = 91;
		if (my_rank == 2 * sub_group)
		{

			MPI_Recv(&tempx[0], 1, MPI_DOUBLE, sub_group, 2 * sub_group, MPI_COMM_WORLD, &status);
			tt[k++] = 92;
			MPI_Recv(&tempy[0], 1, MPI_DOUBLE, sub_group, 2 * sub_group, MPI_COMM_WORLD, &status);
			tt[k++] = 93;

			MPI_Send(&pointx[0], 1, MPI_DOUBLE, 3 * sub_group, 3 * sub_group, MPI_COMM_WORLD);
			tt[k++] = 94;
			MPI_Send(&pointy[0], 1, MPI_DOUBLE, 3 * sub_group, 3 * sub_group, MPI_COMM_WORLD);
			tt[k++] = 95;
		}
		tt[k++] = 96;
		if (my_rank == 3 * sub_group)
		{

			MPI_Recv(&tempx[0], 1, MPI_DOUBLE, 2 * sub_group, 3 * sub_group, MPI_COMM_WORLD, &status);
			tt[k++] = 97;
			MPI_Recv(&tempy[0], 1, MPI_DOUBLE, 2 * sub_group, 3 * sub_group, MPI_COMM_WORLD, &status);
			tt[k++] = 98;

			MPI_Recv(&tempx[1], 1, MPI_DOUBLE, 0, 3 * sub_group, MPI_COMM_WORLD, &status);
			tt[k++] = 99;
			MPI_Recv(&tempy[1], 1, MPI_DOUBLE, 0, 3 * sub_group, MPI_COMM_WORLD, &status);
			tt[k++] = 100;
		}

		MPI_Barrier(MPI_COMM_WORLD);
		tt[k++] = 101;

		
		tt[k++] = 102;
		if (my_rank == sub_group)
		{
			output(my_rank);
			tt[k++] = 103;
		}

		MPI_Barrier(MPI_COMM_WORLD);
		tt[k++] = 104;

		
		tt[k++] = 105;
		if (my_rank == 2 * sub_group)
		{
			output(my_rank);
			tt[k++] = 106;
		}
			
		MPI_Barrier(MPI_COMM_WORLD);
		tt[k++] = 107;

		
		tt[k++] = 108;
		if (my_rank == 3 * sub_group)
		{
			output(my_rank);
			tt[k++] = 109;
		}

		MPI_Barrier(MPI_COMM_WORLD);
		tt[k++] = 110;
		
		if (my_rank == sub_group || my_rank == 2 * sub_group || my_rank == 3 * sub_group)
		{
			MPI_Send(tt, CN, MPI_INT, 0, 6, MPI_COMM_WORLD);
		}
		if (my_rank == 0)
		{
			for (i = 0; i < numpath; i++)
			{
				if (i == 0)
				{
					for (j = 0; j < CN; j++)
						BL[i][j] = tt[j];
				}
				else
				{
					MPI_Recv(BL[i], CN, MPI_INT, i, 6, MPI_COMM_WORLD, &status);
				}
			}
			MPI_Send(BL, numpath*CN, MPI_INT, MASTER, iter + 1, icomm1);
			MPI_Recv(&flag, 1, MPI_INT, MASTER, iter + 2, icomm1, &status);
		}
		MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (flag == 1)
			break;
	}
	
	if (icomm1 != MPI_COMM_NULL) {
		MPI_Barrier(MPI_COMM_WORLD);       
		MPI_Barrier(icomm1);               
		MPI_Comm_disconnect(&icomm1);      
		icomm1 = MPI_COMM_NULL;            
	}

	
	MPI_Finalize();
	return 0;
}
