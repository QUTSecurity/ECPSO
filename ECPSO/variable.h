#include <vector>
#include <algorithm>
#include <map>
#ifndef __VARIABLE_H__
#define __VARIABLE_H__
const int MAXGEN = 100;    // control generations
const int popsize = 30;  // control population size
const int vardim = 6;    // input data dimension

const int numpath = 100;   // number of target paths
const int numprocs = 4;   // number of program processes
const int CN = 300;      // cover node
const int sampsize = 10;  // sampling size
const int selsize = 15;   // selection size

const int MAX_STORED_PARTICLES_PER_PATH = 400; // recommended value, adjustable based on memory

// ... (other extern variables) ...
extern double path_similarity_matrix[numpath][numpath]; // global path similarity matrix
// ... (other constants, structs, etc.) ...
extern int MB[numpath][numprocs][CN]; // declare MB as global variable


const double loser = 0;
const double lprop = 0;
extern int xmax;         // control init pos range
extern int xmin;
extern int vmax;         // control init speed range
extern int vmin;
extern double wmax;      // time-varying weight range
extern double wmin;
extern double c1;        // acceleration coefficient
extern double c2;        // learning factor
extern double thres;     // threshold for non-zero count, keep small
extern double α;
extern double β;
extern double Olim;
extern int NPsQ;
extern int HowUpN;
 

const int  clustersnum = 1000;
const int MASTER = 0;

struct individual
{
    double p[vardim];     // evolved particle
    double speed[vardim];
    double fitness;     // fitness value
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



// Struct: store details of tested particles
struct TestedParticleInfo {
    individual p_info;
    int BL[numprocs][CN]; 
    double fitness;

    TestedParticleInfo() : fitness(0.0) {
        // optional: init BL to 0
        for (int i = 0; i < numprocs; ++i) {
            for (int j = 0; j < CN; ++j) {
                BL[i][j] = 0;
            }
        }
    }

    // constructor: accepts int(*bl)[CN]
    TestedParticleInfo(const individual& particle, int(*bl)[CN], double fit)
        : p_info(particle), fitness(fit) {
        // ======== MODIFICATION START ========
        // direct copy
        for (int i = 0; i < numprocs; ++i) {
            for (int j = 0; j < CN; ++j) {
                BL[i][j] = bl[i][j];
            }
        }
        // ======== MODIFICATION END ========
    }

    // compare operator for descending sort
    bool operator>(const TestedParticleInfo& other) const {
        return fitness > other.fitness;
    }
};





// --- new block: store particle info classified by path fitness ---
struct InitialTestedParticleInfo {
    individual p_info; // particle info
    int bl_array[numprocs][CN]; // path traversal info
    double path_fitness[numpath]; // path fitness [path ID]

    // default constructor
    InitialTestedParticleInfo() {
        // default init for p_info
        // init bl_array to 0
        for (int i = 0; i < numprocs; ++i) {
            for (int j = 0; j < CN; ++j) {
                bl_array[i][j] = 0;
            }
        }
        // init path_fitness to -1
        for (int p = 0; p < numpath; ++p) {
            path_fitness[p] = -1.0;
        }
    }

    // constructor
    InitialTestedParticleInfo(const individual& particle, int(*bl_ptr_to_2d_array)[CN], const double* path_fits)
        : p_info(particle) {
        // copy BL array
        for (int i = 0; i < numprocs; ++i) {
            for (int j = 0; j < CN; ++j) {
                bl_array[i][j] = bl_ptr_to_2d_array[i][j];
            }
        }
        // copy path_fitness array
        for (int p = 0; p < numpath; ++p) {
            path_fitness[p] = path_fits[p];
        }
    }
};

// global variable to store tested particles and their fitness
extern std::vector<InitialTestedParticleInfo> initial_tested_particles_with_fitness;


extern std::vector<TestedParticleInfo> global_tested_particles;
// global map for tested particles per path
extern std::map<int, std::vector<TestedParticleInfo>> path_tested_particles;
// ... (other extern variables) ...
extern double path_similarity_matrix[numpath][numpath]; // global path similarity matrix
// ... (other constants, structs, etc.) ...


#endif // VARIABLE_H_INCLUDED

/* error if redefining */
/* replace with define or const */
// var must be defined in state.cpp when using extern
// state.cpp:
// int a=0; define var and init
// state.h
// extern int a; decl




