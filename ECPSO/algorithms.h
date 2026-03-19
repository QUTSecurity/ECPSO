#ifndef __ALGORITHMS_H__
#define __ALGORITHMS_H__
 
void proposedalgorithm(struct individual populationALL[], double sampling[sampsize + selsize*MAXGEN][vardim], double outputing[sampsize + selsize*MAXGEN], struct  individual population[], double sampling2[sampsize + selsize*MAXGEN][vardim], double outputing2[sampsize + selsize*MAXGEN], struct  individual population2[], int update, MPI_Comm icomm2, double a);

#endif