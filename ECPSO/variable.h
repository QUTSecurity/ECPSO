#include <vector>
#include <algorithm>
#include <map>
#ifndef __VARIABLE_H__
#define __VARIABLE_H__
const int MAXGEN = 100;    
const int popsize = 30;  
const int vardim = 6;    

const int numpath = 100;   
const int numprocs = 4;   
const int CN = 300;      
const int sampsize = 10;  
const int selsize = 15;   

const int MAX_STORED_PARTICLES_PER_PATH = 400; 


extern double path_similarity_matrix[numpath][numpath]; 

extern int MB[numpath][numprocs][CN]; 


const double loser = 0;
const double lprop = 0;
extern int xmax;         
extern int xmin;
extern int vmax;         
extern int vmin;
extern double wmax;      
extern double wmin;
extern double c1;        
extern double c2;        
extern double thres;     
extern double α;
extern double β;
extern double Olim;
extern int NPsQ;
extern int HowUpN;
 

const int  clustersnum = 1000;
const int MASTER = 0;

struct individual
{
    double p[vardim];     
    double speed[vardim];
    double fitness;     
};





template<typename T> std::vector<int> argsort(const std::vector<T>& array)
{
    const int array_len(array.size());
    std::vector<int> array_index(array_len, 0);
    for (int i = 0; i < array_len; ++i)
        array_index[i] = i;

    std::sort(array_index.begin(), array_index.end(),
        [&array](int pos1, int pos2) {return (array[pos1] < array[pos2]); });

    return array_index;
}




struct TestedParticleInfo {
    individual p_info;
    int BL[numprocs][CN]; 
    double fitness;

    TestedParticleInfo() : fitness(0.0) {
        
        for (int i = 0; i < numprocs; ++i) {
            for (int j = 0; j < CN; ++j) {
                BL[i][j] = 0;
            }
        }
    }

    
    TestedParticleInfo(const individual& particle, int(*bl)[CN], double fit)
        : p_info(particle), fitness(fit) {
        
        
        for (int i = 0; i < numprocs; ++i) {
            for (int j = 0; j < CN; ++j) {
                BL[i][j] = bl[i][j];
            }
        }
        
    }

    
    bool operator>(const TestedParticleInfo& other) const {
        return fitness > other.fitness;
    }
};






struct InitialTestedParticleInfo {
    individual p_info; 
    int bl_array[numprocs][CN]; 
    double path_fitness[numpath]; 

    
    InitialTestedParticleInfo() {
        
        
        for (int i = 0; i < numprocs; ++i) {
            for (int j = 0; j < CN; ++j) {
                bl_array[i][j] = 0;
            }
        }
        
        for (int p = 0; p < numpath; ++p) {
            path_fitness[p] = -1.0;
        }
    }

    
    InitialTestedParticleInfo(const individual& particle, int(*bl_ptr_to_2d_array)[CN], const double* path_fits)
        : p_info(particle) {
        
        for (int i = 0; i < numprocs; ++i) {
            for (int j = 0; j < CN; ++j) {
                bl_array[i][j] = bl_ptr_to_2d_array[i][j];
            }
        }
        
        for (int p = 0; p < numpath; ++p) {
            path_fitness[p] = path_fits[p];
        }
    }
};


extern std::vector<InitialTestedParticleInfo> initial_tested_particles_with_fitness;


extern std::vector<TestedParticleInfo> global_tested_particles;

extern std::map<int, std::vector<TestedParticleInfo>> path_tested_particles;

extern double path_similarity_matrix[numpath][numpath]; 



#endif 












