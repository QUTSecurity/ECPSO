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
int tt[CN] = { 0 };//Temporary storage of other segment paths
int BL[numpath][CN] = { 0 };
int k = 0;

/*获得y坐标最小的点,坐标值保存到保存到tempx,tempy数组中*/
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
			/*若y坐标相同,则取x坐标小者。这是为了保证(tempx[0],tempy[0])->(tempx[1],tempy[1])为逆时针方向*/
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

	/*如果y坐标最小的极点只有1点,就复制该点*/
	tempx[1] = tempx[0];
	tempy[1] = tempy[0];
	tt[k++] = 118;

	/*考察极点是否多于1点*/
	for (i = 0; i<n_all; i++)
	{
		tt[k++] = 119;
		tt[k++] = 120;
		if ((pointy[i] == tempy[1]) && (pointx[i]>tempx[1]))
		{
			/*极点多于1点,取x坐标大者,保存最后一个极点的坐标*/
			tempx[1] = pointx[i];
			tempy[1] = pointy[i];
			tt[k++] = 121;
		}
	}
	return;
}


/*获得y坐标最大的点,保存到tempx,tempy数组中*/
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
			/*若y坐标相同,则取x坐标大者。这是为了保证(tempx[0],tempy[0])->(tempx[1],tempy[1])为逆时针方向*/
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
	/*保存YMAX极点坐标*/
	tempx[0] = pointx[index];
	tempy[0] = pointy[index];

	/*如果y坐标最大的极点只有1点,就复制该点*/
	tempx[1] = tempx[0];
	tempy[1] = tempy[0];
	tt[k++] = 129;

	/*考察极点是否多于1点*/
	for (i = 0; i<n_all; i++)
	{
		tt[k++] = 130;
		tt[k++] = 131;
		if ((pointy[i] == tempy[1]) && (pointx[i]<tempx[1]))
		{
			/*极点多于1点,取x坐标小者,保存最后一个极点的坐标*/
			tempx[1] = pointx[i];
			tempy[1] = pointy[i];
			tt[k++] = 132;
		}
	}

	return;
}


/*获得x坐标最小的点,坐标保存到tempx,tempy数组中*/
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
			/*若x坐标相同,则取y坐标大者。这是为了保证(tempx[0],tempy[0])->(tempx[1],tempy[1])为逆时针方向*/
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
	/*保存XMIN极点坐标*/
	tempx[0] = pointx[index];
	tempy[0] = pointy[index];

	/*如果x坐标最小的极点只有1点,就复制该点*/
	tempx[1] = tempx[0];
	tempy[1] = tempy[0];
	tt[k++] = 140;

	/*考察极点是否多于1点*/
	for (i = 0; i<n_all; i++)
	{
		tt[k++] = 141;
		tt[k++] = 142;
		if ((pointx[i] == tempx[1]) && (pointy[i]<tempy[1]))
		{
			/*极点多于1点,取y坐标小者,保存最后一个极点的坐标*/
			tempx[1] = pointx[i];
			tempy[1] = pointy[i];
			tt[k++] = 143;
		}
	}
	return;
}


/*获得x坐标最大的点,其坐标值保存保存到tempx,tempy数组中*/
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
			/*若x坐标相同,则取y坐标小者。这是为了保证(tempx[0],tempy[0])->(tempx[1],tempy[1])为逆时针方向*/
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
	/*保存XMAX极点坐标*/
	tempx[0] = pointx[index];
	tempy[0] = pointy[index];

	/*如果x坐标最大的极点只有1点,就复制该点*/
	tempx[1] = tempx[0];
	tempy[1] = tempy[0];
	tt[k++] = 151;

	/*考察极点是否多于1点*/
	for (i = 0; i<n_all; i++)
	{
		tt[k++] = 152;
		tt[k++] = 153;
		if ((pointx[i] == tempx[1]) && (pointy[i]>tempy[1]))
		{
			/*极点多于1点,取y坐标大者,保存第2个极点的坐标*/
			tempx[1] = pointx[i];
			tempy[1] = pointy[i];
			tt[k++] = 154;
		}
	}
	return;
}


/*temp[x],temp[y]设置两个值的作用是为了有多个最值时，保存首尾两个，后面划分
区域时分别以这两个为端点，中间就没有其他的点了*/

/*确定是属于哪一部分*/
void getincludedvertex(int tag)
{
	double linea, lineb, dist;
	int i, count = 0;
	tt[k++] = 155;

	/*若该边已退化成一个点, 则除此点外无其它点在处理区域内*/
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
		/*根据两端点坐标确定直线的斜率和参数*/
		linea = (tempy[1] - tempy[0]) / (tempx[1] - tempx[0]);
		lineb = tempy[1] - linea*tempx[1];
		tt[k++] = 158;
		/*根据点在直线的哪一侧来确定该顶点是否在处理区域内*/
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
		/*调整各点的位置以加入两个端点*/
		for (i = 0; i<count; i++)
		{
			tt[k++] = 162;
			pointx[count - i] = pointx[count - 1 - i];
			pointy[count - i] = pointy[count - 1 - i];
			tt[k++] = 163;
		}
		/*将两个端点放在两端*/
		pointx[0] = tempx[0];
		pointy[0] = tempy[0];
		pointx[count + 1] = tempx[1];
		pointy[count + 1] = tempy[1];
		count += 2;
		/*更新n_all为要处理的新的表列中的点的数量*/
		n_all = count;
		tt[k++] = 164;
	}
}


/*确定此点的后一点的序号，并放入location[]中*/
void nextindex(int i, int x, int y)
{
	double x1, y1, temp, valuemax = 0;
	int j;
	tt[k++] = 165;

	x1 = pointx[i];
	y1 = pointy[i];

	/*遍历各顶, 将有最小极角点的序号作为自己的Nextindex值*/
	for (j = 0; j<n_all; j++)
	{
		tt[k++] = 166;
		tt[k++] = 167;
		if (((pointx[j] - x1)*x + (pointy[j] - y1)*y)>0)
		{
			/*由极角的余弦或其余角的正弦来判断极角大小*/
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


/*输出点的序列*/
void output(int rank)
{
	int j = 0, index, temp;
	double xtem, ytem;
	int flag = 0;
	tt[k++] = 173;

	//printf("输出的是第%d部分的点\n", rank / sub_group);
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

	/*将nextindex反向,使输出的点按逆时针方向,以保证各行主处理器输出顺序相邻*/
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

	/*按照索引顺序输出各点*/
	j = n_all - 1;
	while (location[j] != -1)
	{
		tt[k++] = 183;
		tt[k++] = 184;
		if (!((j == n_all - 1) && (flag == 0)))
		{
			//printf("%.2lf,%.2lf\n", pointx[j], pointy[j]);
			xtem = pointx[j];
			ytem = pointy[j];
			tt[k++] = 185;
		}
		j = location[j];
		tt[k++] = 186;
	}

	/*判断是不是重合*/
	tt[k++] = 187;
	if (!((j == n_all - 1) && (flag == 0)))
	{
		//printf("%.2lf,%.2lf\n", pointx[j], pointy[j]);
		xtem = pointx[j];
		ytem = pointy[j];
		tt[k++] = 188;
	}

	tt[k++] = 189;
	if ((rank == 3 * sub_group) && !((xtem == tempx[1]) && (ytem == tempy[1])))
		tt[k++] = 190;
	//printf("%.2lf,%.2lf\n", tempx[1], tempy[1]);
	return;
}


int main(int argc, char *argv[])
{
	int i, j, flag = 0;
	int n;             //定义变量防溢出
	int looptimes = 0, iter = 0;
	int points[vardim] = { 0 };
	MPI_Comm icomm1;
	MPI_Init(&argc, &argv);
	MPI_Comm_get_parent(&icomm1);
	MPI_Comm_size(MPI_COMM_WORLD, &group_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	/*
	//本程序至少要4个处理器才能正常执行
	if (group_size < 4)
	{
		if (my_rank == 0)
		{
			printf("Need 4 or more processors to run!\n");
		}
		MPI_Finalize();
		exit(0);
	}

	if (my_rank == 0)
	{
		printf("please input all the vertexes!\nfirst is number!\n");
		printf("please input the Number:");
		scanf("%d", &n_all);
		printf("please input the vertex:\n");
		for (i = 0; i<n_all; i++)
		{
			scanf("%lf", &pointx[i]);
			scanf("%lf", &pointy[i]);
		}
	}
	//处理器0接收输入的顶点坐标
	*/
	sub_group = group_size / 4;
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	if (my_rank == 0)
	{
		MPI_Recv(&looptimes, 1, MPI_INT, MASTER, 11, icomm1, &status);
	}
	MPI_Bcast(&looptimes, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	//cout << "looptimes=" << looptimes << endl;
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
		///////////////////////////////////////////////////////////////////////////////////////////////////
		tt[k++] = 1;
		MPI_Bcast(&n_all, 1, MPI_INT, 0, MPI_COMM_WORLD);
		tt[k++] = 2;
		/*广播顶点总数*/
		/*第1行的行主处理器,向第2,3,4行的行主处理器发送顶点坐标，对应于算法17.6步骤(1.1)*/
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

		/*计算极点,并把他们的坐标分别存储在4个行主处理器中*/
		/*第1行的行主处理器计算YMIN*/
		tt[k++] = 11;
		if (my_rank == 0)
		{
			getymin();
			tt[k++] = 12;
		}
		/*第2行的行主处理器计算XMAX*/
		tt[k++] = 13;
		if (my_rank == sub_group)
		{
			getxmax();
			tt[k++] = 14;
		}
		/*第3行的行主处理器计算YMAX*/
		tt[k++] = 15;
		if (my_rank == 2 * sub_group)
		{
			getymax();
			tt[k++] = 16;
		}
		/*第4行的行主处理器计算XMIN*/
		tt[k++] = 17;
		if (my_rank == 3 * sub_group)
		{
			getxmin();
			tt[k++] = 18;
		}

		MPI_Barrier(MPI_COMM_WORLD);
		tt[k++] = 19;

		/*将四条由极点组成的边存储到每一行的行主处理器上*/
		tt[k++] = 20;
		if ((my_rank>0) && (my_rank<4 * sub_group) && (my_rank%sub_group == 0))
		{
			/*各行主处理器发送消息给进程0, 告知极点坐标*/
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
				/*进程0将相关的极点坐标发送给各行主处理器*/
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
				/*各行主处理器接收进程0的消息,得到并存储由极点构成的边*/
				MPI_Recv(&tempx[1], 1, MPI_DOUBLE, 0, my_rank, MPI_COMM_WORLD, &status);
				tt[k++] = 33;
				MPI_Recv(&tempy[1], 1, MPI_DOUBLE, 0, my_rank, MPI_COMM_WORLD, &status);
				tt[k++] = 34;
			}
		}

		MPI_Barrier(MPI_COMM_WORLD);
		tt[k++] = 35;

		/*确定四边形中的顶点，并将其余顶点归入四个三角区中*/
		/*四个行主处理器同时判断顶点是否处于自身所在的区域 */

		/*保留属于该区域内的点，每个行主处理器得到要处理的新的表列*/
		/*下面得到属于此范围的点*/
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
				/*各行主处理器向该行的其他处理器发送新的表列中的顶点数量*/
				MPI_Send(&n_all, 1, MPI_INT, my_rank + i, i, MPI_COMM_WORLD);
				tt[k++] = 47;
			}
		}
		else
		{
			tt[k++] = 48;
			if (my_rank < 4 * sub_group)
			{
				/*各处理器从该行的主处理器接收新的表列中的顶点数量*/
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
				/*各行主处理器向该行的其他处理器发送新的表列中的顶点坐标*/
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

		/*每一行上的处理器对属于同一区域的点计算极角；将有最小极角点的序号作为自己的Nexindex值建立各进程中点的nextindex索引*/
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

		/*将个进程中点的nextindex索引发送到中心结点*/
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

		/*每个行主处理器按照点的Nextindex索引输出自己处理器上的点)*/
		tt[k++] = 80;
		if (my_rank == 0)
		{
			/*第1行输出各顶点的坐标*/
			output(my_rank);
			tt[k++] = 81;
			/*将最后的一个点传送给第2行主处理器判断是不是相同*/
			MPI_Send(&pointx[0], 1, MPI_DOUBLE, sub_group, sub_group, MPI_COMM_WORLD);
			tt[k++] = 82;
			MPI_Send(&pointy[0], 1, MPI_DOUBLE, sub_group, sub_group, MPI_COMM_WORLD);
			tt[k++] = 83;
			/*将起始点传送给第4行主处理器判断是不是相同*/
			MPI_Send(&pointx[n_all - 1], 1, MPI_DOUBLE, 3 * sub_group, 3 * sub_group, MPI_COMM_WORLD);
			tt[k++] = 84;
			MPI_Send(&pointy[n_all - 1], 1, MPI_DOUBLE, 3 * sub_group, 3 * sub_group, MPI_COMM_WORLD);
			tt[k++] = 85;
		}
		tt[k++] = 86;
		if (my_rank == sub_group)
		{
			/*接收第1行主处理器发送过来的最后一个点,也即本行起始点,以判断是不是相同*/
			MPI_Recv(&tempx[0], 1, MPI_DOUBLE, 0, sub_group, MPI_COMM_WORLD, &status);
			tt[k++] = 87;
			MPI_Recv(&tempy[0], 1, MPI_DOUBLE, 0, sub_group, MPI_COMM_WORLD, &status);
			tt[k++] = 88;
			/*将最后的一个点传送给第3行主处理器判断是不是相同*/
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

		/*第2行输出各顶点的坐标*/
		tt[k++] = 102;
		if (my_rank == sub_group)
		{
			output(my_rank);
			tt[k++] = 103;
		}

		MPI_Barrier(MPI_COMM_WORLD);
		tt[k++] = 104;

		/*第3行输出各顶点的坐标*/
		tt[k++] = 105;
		if (my_rank == 2 * sub_group)
		{
			output(my_rank);
			tt[k++] = 106;
		}
			
		MPI_Barrier(MPI_COMM_WORLD);
		tt[k++] = 107;

		/*第4行输出各顶点的坐标*/
		tt[k++] = 108;
		if (my_rank == 3 * sub_group)
		{
			output(my_rank);
			tt[k++] = 109;
		}

		MPI_Barrier(MPI_COMM_WORLD);
		tt[k++] = 110;
		///////////////////////////////////////////////////////////////////////////////////
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
	// 在主程序结束时，断开通信器
	if (icomm1 != MPI_COMM_NULL) {
		MPI_Barrier(MPI_COMM_WORLD);       // 进行全局同步
		MPI_Barrier(icomm1);               // 父进程同步
		MPI_Comm_disconnect(&icomm1);      // 断开主通信器
		icomm1 = MPI_COMM_NULL;            // 清空句柄
	}

	/*
	if (my_rank == 0)
	{
		cout << "BL=" << endl;
		cout << "{";
		for (i = 0; i < numpath; i++)
		{
			cout << "{";
			for (j = 0; j < CN; j++)
			{
				cout << BL[i][j] << ",";
				if (j % (CN / 8) == 0 && j != 0)
				{
					cout << "\n";
				}
			}
			cout << "}";
			if (i < numpath - 1)
			{
				cout << ",";
				cout << "\n";
			}
		}
		cout << "}" << endl;
	}
	*/
	MPI_Finalize();
	return 0;
}
