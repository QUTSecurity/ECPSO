#define MPICH_SKIP_MPICXX
#include <mpi.h>
#include "variable.h"
#include "fitness.h"
#include "latin.h"
#include "algorithms.h"
#include <iostream>
#include <cstdio> 
#include <cstdlib>
#include <cmath>

#include <cstring>
#include <ctime>
#include <vector>
#include <algorithm>
#include <numeric>
#include <map>
#include <tuple>

using namespace std;






enum class GroupingMode {
	Auto,
	StaticOnly,
	DynamicOnly
};

const GroupingMode GROUPING_MODE = GroupingMode::Auto;
const double THETA_SWITCH = 0.8;
const double THETA_GROUP = 0.6;
static const char* groupingModeName(GroupingMode mode) {
	switch (mode) {
	case GroupingMode::Auto: return "auto";
	case GroupingMode::StaticOnly: return "static_only";
	case GroupingMode::DynamicOnly: return "dynamic_only";
	default: return "unknown";
	}
}
static void printSimilarityStats(const char* name, double M[][numpath], double theta) {
	double max_sim = 0.0;
	double avg_sim = 0.0;
	int above_theta = 0;
	int cnt = 0;
	for (int i = 0; i < numpath; ++i) {
		for (int j = i + 1; j < numpath; ++j) {
			double v = M[i][j];
			avg_sim += v;
			if (v > max_sim) max_sim = v;
			if (v >= theta) above_theta++;
			cnt++;
		}
	}
	if (cnt > 0) avg_sim /= cnt;
	printf("[GROUPING] %s stats: max=%.4f, avg=%.4f, pairs>=theta=%d/%d (theta=%.2f)\n",
		name, max_sim, avg_sim, above_theta, cnt, theta);
	fflush(stdout);
}
static void computeStaticSimilarityMatrix(double M[][numpath]) {
	for (int i = 0; i < numpath; ++i) {
		M[i][i] = 1.0;
		for (int j = i + 1; j < numpath; ++j) {
			double s = static_cast<double>(Sim(MB[i], MB[j]));
			M[i][j] = s;
			M[j][i] = s;
		}
	}
}
static bool allStaticBelowThreshold(double M[][numpath], double theta) {
	for (int i = 0; i < numpath; ++i) {
		for (int j = i + 1; j < numpath; ++j) {
			if (M[i][j] >= theta) return false;
		}
	}
	return true;
}
static double computeSpearmanCorrelation(const double* f_i, const double* f_j, int n) {
	
	auto rank_average = [&](const double* f, vector<double>& r) {
		vector<int> idx(n);
		iota(idx.begin(), idx.end(), 0);
		stable_sort(idx.begin(), idx.end(), [&](int a, int b) { return f[a] < f[b]; });
		int k = 0;
		while (k < n) {
			int start = k;
			double val = f[idx[k]];
			while (k + 1 < n && f[idx[k + 1]] == val) {
				k++;
			}
			int end = k;
			double avg_rank = 0.5 * (start + end);
			for (int t = start; t <= end; ++t) {
				r[idx[t]] = avg_rank;
			}
			k++;
		}
		};
	vector<double> r_i(n), r_j(n);
	rank_average(f_i, r_i);
	rank_average(f_j, r_j);
	double mean_ri = 0.0, mean_rj = 0.0;
	for (int k = 0; k < n; ++k) {
		mean_ri += r_i[k];
		mean_rj += r_j[k];
	}
	mean_ri /= n;
	mean_rj /= n;
	double cov_val = 0.0, var_ri = 0.0, var_rj = 0.0;
	for (int k = 0; k < n; ++k) {
		double di = r_i[k] - mean_ri;
		double dj = r_j[k] - mean_rj;
		cov_val += di * dj;
		var_ri += di * di;
		var_rj += dj * dj;
	}
	double denom = sqrt(var_ri * var_rj);
	return (denom < 1e-10) ? 0.0 : cov_val / denom;
}
static void computeDynamicSimilarityMatrix(double Md[][numpath], double outputs[][sampsize]) {
	for (int i = 0; i < numpath; ++i) {
		Md[i][i] = 1.0;
		for (int j = i + 1; j < numpath; ++j) {
			double rho = computeSpearmanCorrelation(outputs[i], outputs[j], sampsize);
			double s_dyn = max(0.0, (rho + 1.0) / 2.0);
			Md[i][j] = s_dyn;
			Md[j][i] = s_dyn;
		}
	}
}
static vector<vector<int>> greedyGrouping(double M_final[][numpath], double theta) {
	vector<vector<int>> groups;
	vector<bool> assigned(numpath, false);
	vector<tuple<double, int, int>> pairs;
	for (int i = 0; i < numpath; ++i) {
		for (int j = i + 1; j < numpath; ++j) {
			pairs.emplace_back(M_final[i][j], i, j);
		}
	}
	sort(pairs.rbegin(), pairs.rend());
	for (const auto& pr : pairs) {
		double sim = get<0>(pr);
		int pi = get<1>(pr);
		int pj = get<2>(pr);
		if (sim < theta) break;
		if (!assigned[pi] && !assigned[pj]) {
			groups.push_back({ pi, pj });
			assigned[pi] = true;
			assigned[pj] = true;
		}
		else if (assigned[pi] != assigned[pj]) {
			int assigned_one = assigned[pi] ? pi : pj;
			int u = assigned[pi] ? pj : pi;
			for (auto& g : groups) {
				if (find(g.begin(), g.end(), assigned_one) != g.end()) {
					double avgSim = 0.0;
					for (int v : g) avgSim += M_final[u][v];
					avgSim /= g.size();
					if (avgSim >= theta) {
						g.push_back(u);
						assigned[u] = true;
					}
					break;
				}
			}
		}
	}
	for (int i = 0; i < numpath; ++i) {
		if (!assigned[i]) groups.push_back({ i });
	}
	return groups;
}
int main(int argc, char** argv)
{
	
	int path_start_index = 0;
	int path_end_index = numpath - 1;
	if (argc == 3) {
		path_start_index = atoi(argv[1]);
		path_end_index = atoi(argv[2]);
		cout << "[INFO] Batch mode: Processing paths [" << path_start_index << ", " << path_end_index << "]" << endl;
	}
	else {
		cout << "[INFO] Full mode: Processing all " << numpath << " paths" << endl;
		cout << "[INFO] Usage for batch mode: .\\FERPSO_Ctest.exe <start_path_index> <end_path_index>" << endl;
	}
	if (path_start_index < 0 || path_end_index >= numpath || path_start_index > path_end_index) {
		cerr << "[ERROR] Invalid path range. Exiting." << endl;
		return 1;
	}
	int b = 0, i = 0, j = 0, k = 0, l = 0, t = 0, loserPop = 0;
	int iter = 0, count[numpath] = { 0 };
	int update = 1;
	int flag[numpath] = { 0 };
	int testinput[vardim] = { 1, 1, 1 };
	double r1, r2, w;
	int tempint = 0;
	double temp = 0.0, temp1 = 0.0, temp2 = 0.0, temp3 = 0.0, temp4 = 0.0;
	double perbest[popsize][vardim] = { 0 };
	double perbestfit[popsize] = { 0 };
	double globest[vardim] = { 0 };
	double globestfit = 0.0;
	double gloworst[vardim] = { 0 };
	double gloworstfit = 10000.0;
	int numRBFN = 2;
	double a = 0.5;
	double Fgood1 = 0.0;
	double Fgood2 = 0.0;
	int numCountG1 = 0;
	int numCountG2 = 0;
	struct individual personLeft;
	struct individual personRight;
	int maxUp = 500, minUp = 50;
	int(*BL)[numprocs][CN] = new int[popsize][numprocs][CN];
	memset(BL, 0, sizeof(int) * popsize * numprocs * CN);
	double(*samples)[vardim] = new double[sampsize][vardim];
	memset(samples, 0, sizeof(double) * sampsize * vardim);
	double(*outputs)[sampsize] = new double[numpath][sampsize];
	memset(outputs, 0, sizeof(double) * numpath * sampsize);
	double(*samples2)[vardim] = new double[sampsize][vardim];
	memset(samples2, 0, sizeof(double) * sampsize * vardim);
	double(*outputs2)[sampsize] = new double[numpath][sampsize];
	memset(outputs2, 0, sizeof(double) * numpath * sampsize);
	int(*sampBL)[numprocs][CN] = new int[sampsize][numprocs][CN];
	memset(sampBL, 0, sizeof(int) * sampsize * numprocs * CN);
	double(*sampling)[sampsize + selsize * MAXGEN][vardim] = new double[numpath][sampsize + selsize * MAXGEN][vardim];
	memset(sampling, 0, sizeof(double) * numpath * (sampsize + selsize * MAXGEN) * vardim);
	double(*outputing)[sampsize + selsize * MAXGEN] = new double[numpath][sampsize + selsize * MAXGEN];
	memset(outputing, 0, sizeof(double) * numpath * (sampsize + selsize * MAXGEN));
	double(*sampling2)[sampsize + selsize * MAXGEN][vardim] = new double[numpath][sampsize + selsize * MAXGEN][vardim];
	memset(sampling2, 0, sizeof(double) * numpath * (sampsize + selsize * MAXGEN) * vardim);
	double(*outputing2)[sampsize + selsize * MAXGEN] = new double[numpath][sampsize + selsize * MAXGEN];
	memset(outputing2, 0, sizeof(double) * numpath * (sampsize + selsize * MAXGEN));
	vector<int> index;
	vector<int> index1;
	vector<int> index2;
	vector<int> index3;
	vector<double> estfit;
	vector<double> estfit1;
	vector<double> estfit2;
	vector<double> estfit3;
	double starttime, endtime, consumedtimes[numpath] = { 0 };
	struct individual population[popsize];
	struct individual population2[popsize];
	struct individual populationALL[2 * popsize];
	struct individual superior_inds[selsize];
	int testdata[numpath][vardim] = { 0 };
	
	int program_executions = 0;
	
	
	
	const size_t EXPPOOL_MAX_SIZE = MAX_STORED_PARTICLES_PER_PATH * numpath; 
	std::vector<TestedParticleInfo> global_particles;
	global_particles.reserve(EXPPOOL_MAX_SIZE); 
	size_t pool_write_idx = 0;  
	bool pool_full = false;     
	int cu_sampling[numpath] = { 0 }; 
#ifdef _WIN32
	char programpath[100] = "../convex/convex.exe"; 
	char rbfnpath[100] = "../rbfn/rbfn.exe";
#else
	char programpath[100] = "../convex/convex"; 
	char rbfnpath[100] = "../rbfn/rbfn";
#endif
	int myid, size;
	srand((unsigned)time(NULL));
	MPI_Status status;
	MPI_Comm icomm2;
	MPI_Comm icomm1;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	starttime = MPI_Wtime();
	
	latin_hypercube(samples);
	MPI_Comm_spawn(programpath, MPI_ARGV_NULL, numprocs, MPI_INFO_NULL, MASTER, MPI_COMM_SELF, &icomm1, MPI_ERRCODES_IGNORE);
	tempint = sampsize;
	MPI_Send(&tempint, 1, MPI_INT, 0, 11, icomm1);
	for (i = 0; i < sampsize; i++)
	{
		for (j = 0; j < vardim; j++)
		{
			testinput[j] = samples[i][j];
		}
		MPI_Send(testinput, vardim, MPI_INT, 0, i, icomm1);
		program_executions++; 
		MPI_Recv(sampBL[i], numprocs * CN, MPI_INT, 0, i + 1, icomm1, &status);
		for (b = 0; b < numpath; b++)
		{
			outputs[b][i] = Sim(sampBL[i], MB[b]);
		}
		
		{
			individual particle;
			for (int dim = 0; dim < vardim; ++dim) {
				particle.p[dim] = samples[i][dim];
				particle.speed[dim] = 0.0;
			}
			particle.fitness = outputs[0][i];
			int(*bl)[CN] = &sampBL[i][0];
			
			if (!pool_full) {
				global_particles.emplace_back(particle, bl, particle.fitness);
				if (global_particles.size() >= EXPPOOL_MAX_SIZE) pool_full = true;
			}
			else {
				global_particles[pool_write_idx] = TestedParticleInfo(particle, bl, particle.fitness);
			}
			pool_write_idx = (pool_write_idx + 1) % EXPPOOL_MAX_SIZE;
		}
		MPI_Send(&flag[0], 1, MPI_INT, 0, i + 2, icomm1);
	}
	for (b = 0; b < numpath; b++)
	{
		for (i = 0; i < sampsize; i++)
		{
			for (j = 0; j < vardim; j++)
			{
				sampling[b][i][j] = samples[i][j];
				sampling2[b][i][j] = samples[i][j];
			}
			outputing[b][i] = outputs[b][i];
			outputing2[b][i] = outputs[b][i];
		}
	}
	double static_sim[numpath][numpath] = { 0 };
	double dynamic_sim[numpath][numpath] = { 0 };
	double final_sim[numpath][numpath] = { 0 };
	computeStaticSimilarityMatrix(static_sim);
	printf("[GROUPING] mode=%s, theta_switch=%.2f, theta_group=%.2f\n",
		groupingModeName(GROUPING_MODE), THETA_SWITCH, THETA_GROUP);
	printSimilarityStats("static", static_sim, THETA_GROUP);
	bool use_dynamic = false;
	if (GROUPING_MODE == GroupingMode::DynamicOnly) {
		use_dynamic = true;
	}
	else if (GROUPING_MODE == GroupingMode::StaticOnly) {
		use_dynamic = false;
	}
	else {
		
		use_dynamic = allStaticBelowThreshold(static_sim, THETA_SWITCH);
	}
	if (use_dynamic) {
		computeDynamicSimilarityMatrix(dynamic_sim, outputs);
		memcpy(final_sim, dynamic_sim, sizeof(final_sim));
		printSimilarityStats("dynamic", dynamic_sim, THETA_GROUP);
		printf("[GROUPING] strategy=dynamic\n");
	}
	else {
		memcpy(final_sim, static_sim, sizeof(final_sim));
		printf("[GROUPING] strategy=static\n");
	}
	vector<vector<int>> path_groups = greedyGrouping(final_sim, THETA_GROUP);
	printf("[GROUPING] groups=%d\n", (int)path_groups.size());
	for (size_t gi = 0; gi < path_groups.size(); ++gi) {
		printf("[GROUPING] group %d: [", (int)gi);
		for (size_t gj = 0; gj < path_groups[gi].size(); ++gj) {
			printf("%d%s", path_groups[gi][gj], (gj + 1 < path_groups[gi].size()) ? ", " : "");
		}
		printf("]\n");
	}
	fflush(stdout);
	
	
	const int STAGNATION_W = 5;
	const double STAGNATION_EPS = 1e-6;
	const int ESCAPE_K = 15;
	vector<vector<double>> group_best_hist(path_groups.size());
	vector<int> group_escape_stage(path_groups.size(), 0); 
	for (size_t group_idx = 0; group_idx < path_groups.size(); ++group_idx)
	{
		const vector<int>& current_group = path_groups[group_idx];
		int group_anchor = -1;
		for (size_t member_idx = 0; member_idx < current_group.size(); ++member_idx) {
			int candidate = current_group[member_idx];
			if (candidate >= path_start_index && candidate <= path_end_index) {
				group_anchor = candidate;
				break;
			}
		}
		if (group_anchor < 0) {
			continue;
		}
		b = group_anchor;
		bool group_finished = false;
		
		
		
		{
			const double EPSILON_TRANSFER = thres; 
			for (size_t gm = 0; gm < current_group.size(); ++gm) {
				int path_id = current_group[gm];
				if (path_id < path_start_index || path_id > path_end_index || flag[path_id] == 1) {
					continue;
				}
				for (size_t pi = 0; pi < global_particles.size(); ++pi) {
					double sim_val = Sim(global_particles[pi].BL, MB[path_id]);
					if (fabs(sim_val - 1.0) < EPSILON_TRANSFER) {
						
						flag[path_id] = 1;
						consumedtimes[path_id] = MPI_Wtime() - starttime;
						count[path_id] = 0; 
						for (int dim = 0; dim < vardim; ++dim) {
							testdata[path_id][dim] = global_particles[pi].p_info.p[dim];
						}
						int current_solved = 0;
						for (int c = 0; c < numpath; ++c) {
							if (flag[c] == 1) current_solved++;
						}
						printf("[REUSE_SUCCESS] PathID=%d, GroupID=%d, Time=%.4fs, Complete=%.2f%% (%d/%d)\n",
							path_id, (int)group_idx, consumedtimes[path_id],
							(double)current_solved / numpath * 100.0, current_solved, numpath);
						fflush(stdout);
						break; 
					}
				}
			}
			
			bool all_solved_by_transfer = true;
			for (size_t gm = 0; gm < current_group.size(); ++gm) {
				int path_id = current_group[gm];
				if (path_id >= path_start_index && path_id <= path_end_index && flag[path_id] == 0) {
					all_solved_by_transfer = false;
					break;
				}
			}
			if (all_solved_by_transfer) {
				printf("[TRANSFER_SKIP] GroupID=%d all paths covered by pool, skipping search\n", (int)group_idx);
				fflush(stdout);
				continue; 
			}
		}
		
		
		cu_sampling[b] = sampsize; 
		{
			const double FITNESS_THRESHOLD = 0.25;  
			const int    TRANSFER_TOTAL = 0;      
			const int    QUOTA_TOP = 100;            
			const int    TRANSFER_PER_PATH = 40;     
			
			struct TransferCandidate { const TestedParticleInfo* particle; double max_fitness, min_fitness; };
			std::vector<TransferCandidate> candidates;
			candidates.reserve(global_particles.size());
			for (const auto& tp : global_particles) {
				double mx = -1e100, mn = 1e100;
				for (int path_id : current_group) {
					if (flag[path_id] == 1 || path_id < 0 || path_id >= numpath) continue;
					double f = Sim(tp.BL, MB[path_id]);
					if (f > mx) mx = f;
					if (f < mn) mn = f;
				}
				if (mx >= FITNESS_THRESHOLD) candidates.push_back({ &tp, mx, mn });
			}
			if (candidates.empty()) {
				printf("[TRANSFER_SUMMARY] Group %d: screened=0, transferred=0\n", (int)group_idx);
				fflush(stdout);
			}
			else {
				
				std::sort(candidates.begin(), candidates.end(),
					[](const TransferCandidate& a, const TransferCandidate& b) {
						return a.max_fitness > b.max_fitness;
					});
				
				std::map<int, int> path_contrib;
				for (int path_id : current_group) path_contrib[path_id] = 0;
				
				auto choose_owner = [&](const TransferCandidate& c)->int {
					int owner = -1; double best = -1e100;
					for (int path_id : current_group) {
						if (flag[path_id] == 1) continue;
						if (path_id < 0 || path_id >= numpath) continue;
						if (path_contrib[path_id] >= TRANSFER_PER_PATH) continue;
						double f = Sim(c.particle->BL, MB[path_id]);
						if (f > FITNESS_THRESHOLD && f > best) { best = f; owner = path_id; }
					}
					return owner;
					};
				
				auto append = [&](const TransferCandidate& c, int owner)->bool {
					if (owner == -1) return false;
					if (cu_sampling[b] >= (sampsize + selsize * MAXGEN)) return false;
					for (int d = 0; d < vardim; ++d) {
						sampling[b][cu_sampling[b]][d] = c.particle->p_info.p[d];
						sampling2[b][cu_sampling[b]][d] = c.particle->p_info.p[d];
					}
					outputing[b][cu_sampling[b]] = c.min_fitness;  
					outputing2[b][cu_sampling[b]] = c.max_fitness; 
					cu_sampling[b]++; path_contrib[owner]++; return true;
					};
				int transferred = 0;
				
				for (int idx = 0; idx < (int)candidates.size()
					&& transferred < QUOTA_TOP
					&& transferred < TRANSFER_TOTAL; ++idx) {
					int owner = choose_owner(candidates[idx]);
					if (append(candidates[idx], owner)) transferred++;
				}
				
				if (transferred < TRANSFER_TOTAL) {
					int remain_quota = TRANSFER_TOTAL - transferred;
					int step = std::max(1, (int)candidates.size() / std::max(1, remain_quota));
					for (int idx = 0; idx < (int)candidates.size() && transferred < TRANSFER_TOTAL; idx += step) {
						int owner = choose_owner(candidates[idx]);
						if (append(candidates[idx], owner)) transferred++;
					}
				}
				printf("[TRANSFER_SUMMARY] Group %d: candidates=%d, transferred=%d, train_len=%d\n",
					(int)group_idx, (int)candidates.size(), transferred, cu_sampling[b]);
				fflush(stdout);
			}
		}
		for (i = 0; i < popsize; i++)
		{
			for (j = 0; j < vardim; j++)
			{
				population[i].p[j] = fabs((float)(rand() % (xmax - xmin)) + xmin);        
				population[i].speed[j] = fabs((float)(rand() % (vmax - vmin + 1)) + vmin);
				population2[i].p[j] = fabs((float)(rand() % (xmax - xmin)) + xmin);        
				population2[i].speed[j] = fabs((float)(rand() % (vmax - vmin + 1)) + vmin);
			}
		}
		MPI_Comm_spawn(programpath, MPI_ARGV_NULL, numprocs, MPI_INFO_NULL, MASTER, MPI_COMM_SELF, &icomm1, MPI_ERRCODES_IGNORE);
		tempint = selsize * MAXGEN;
		MPI_Send(&tempint, 1, MPI_INT, 0, 11, icomm1);
		MPI_Comm_spawn(rbfnpath, MPI_ARGV_NULL, numRBFN, MPI_INFO_NULL, MASTER, MPI_COMM_SELF, &icomm2, MPI_ERRCODES_IGNORE);
		MPI_Send(&MAXGEN, 1, MPI_INT, 0, 22, icomm2);
		for (iter = 0; iter < MAXGEN; iter++)
		{
			
			a = 1 / (1 + exp(-iter));
			memset(BL, 0, sizeof(int) * popsize * numprocs * CN);
			MPI_Send(&update, 1, MPI_INT, 0, iter, icomm2);
			proposedalgorithm(populationALL, sampling[b], outputing[b], population, sampling2[b], outputing2[b], population2, update, icomm2, a);
			numCountG1 = 0;
			numCountG2 = 0;
			for (i = 0; i < popsize; i++)
			{
				if (population[i].fitness >= Olim) {
					numCountG1++;
				}
				if (population2[i].fitness >= Olim) {
					numCountG2++;
				}
			}
			Fgood1 = numCountG1 / popsize;
			Fgood2 = numCountG2 / popsize;
			if (Fgood1 == 0 && Fgood2 == 0) {
				estfit.clear();
				for (i = 0; i < 2 * popsize; i++)
				{
					estfit.push_back(populationALL[i].fitness);
				}
				index.clear();
				index = argsort(estfit);
				for (i = 0; i < selsize; i++)
				{
					for (j = 0; j < vardim; j++)
					{
						superior_inds[i].p[j] = populationALL[index[abs(2 * popsize - i - 1)]].p[j];
					}
					superior_inds[i].fitness = populationALL[index[abs(2 * popsize - i - 1)]].fitness;
				}
			}
			else {
				int F1 = (Fgood1 / (Fgood1 + Fgood2)) * selsize;
				estfit3.clear();
				for (i = 0; i < popsize; i++)
				{
					estfit3.push_back(population[i].fitness);
				}
				index3.clear();
				index3 = argsort(estfit3);
				for (i = 0; i < F1; i++)
				{
					for (j = 0; j < vardim; j++)
					{
						superior_inds[i] = population[index3[abs(popsize - i - 1)]];
					}
				}
				estfit3.clear();
				for (i = 0; i < popsize; i++)
				{
					estfit3.push_back(population2[i].fitness);
				}
				index3.clear();
				index3 = argsort(estfit3);
				for (i = F1; i < selsize; i++)
				{
					for (j = 0; j < vardim; j++)
					{
						superior_inds[i] = population2[index3[abs(popsize - i - 1)]];
					}
				}
			}
			
			
			
			
			for (i = 0; i < selsize; i++)
			{
				program_executions++; 
				
				for (j = 0; j < vardim; j++)
				{
					testinput[j] = superior_inds[i].p[j];
				}
				MPI_Send(testinput, vardim, MPI_INT, 0, (iter * selsize) + i, icomm1);
				MPI_Recv(BL[i], numprocs * CN, MPI_INT, 0, (iter * selsize) + i + 1, icomm1, &status);
				double predicted_fitness = superior_inds[i].fitness;
				superior_inds[i].fitness = Sim(BL[i], MB[b]);
				double prediction_error = fabs(predicted_fitness - superior_inds[i].fitness);
				(void)prediction_error;
				
				
				{
					if (!pool_full) {
						global_particles.emplace_back(superior_inds[i], BL[i], superior_inds[i].fitness);
						if (global_particles.size() >= EXPPOOL_MAX_SIZE) pool_full = true;
					}
					else {
						global_particles[pool_write_idx] = TestedParticleInfo(superior_inds[i], BL[i], superior_inds[i].fitness);
					}
					pool_write_idx = (pool_write_idx + 1) % EXPPOOL_MAX_SIZE;
				}
				
				int current_solved = 0;
				for (int c = 0; c < numpath; ++c) {
					if (flag[c] == 1) current_solved++;
				}
				double completion_rate = (double)current_solved / numpath * 100.0;
				for (size_t gm = 0; gm < current_group.size(); ++gm) {
					int path_id = current_group[gm];
					if (path_id < path_start_index || path_id > path_end_index || flag[path_id] == 1) {
						continue;
					}
					double hit_fitness = Sim(BL[i], MB[path_id]);
					if (fabs(hit_fitness - 1.0) < thres) {
						flag[path_id] = 1;
						consumedtimes[path_id] = MPI_Wtime() - starttime;
						count[path_id] = iter + 1;
						int solved_after = current_solved + 1;
						double completion_rate_after = (double)solved_after / numpath * 100.0;
						printf("[PSO_SUCCESS] PathID=%d, GroupID=%d, Time=%.4fs, Iter=%d, Complete=%.2f%% (%d/%d)\n",
							path_id, (int)group_idx, consumedtimes[path_id], iter + 1, completion_rate_after, solved_after, numpath);
						fflush(stdout);
						for (j = 0; j < vardim; j++) {
							testdata[path_id][j] = superior_inds[i].p[j];
						}
					}
				}
				int anchor_done = 0;
				MPI_Send(&anchor_done, 1, MPI_INT, 0, (iter * selsize) + i + 2, icomm1);
			}

			
			{
				double current_best = -1e100;
				for (int si = 0; si < selsize; ++si) {
					if (superior_inds[si].fitness > current_best) {
						current_best = superior_inds[si].fitness;
					}
				}

				auto& hist = group_best_hist[group_idx];
				hist.push_back(current_best);
				bool stagnated = false;
				double delta_w = 0.0;
				
				if ((int)hist.size() > STAGNATION_W) {
					delta_w = fabs(hist.back() - hist[(int)hist.size() - 1 - STAGNATION_W]);
					stagnated = (delta_w < STAGNATION_EPS);
				}

				if (stagnated) {
					printf("[STAGNATION_DETECT] group=%d, iter=%d, best=%.6f, delta_w=%.8f\n",
						(int)group_idx, iter, current_best, delta_w);
					fflush(stdout);

					if (group_escape_stage[group_idx] == 0) {
						
						
						vector<int> unsolved_paths;
						for (size_t gm = 0; gm < current_group.size(); ++gm) {
							int path_id = current_group[gm];
							if (path_id >= path_start_index && path_id <= path_end_index && flag[path_id] == 0) {
								unsolved_paths.push_back(path_id);
							}
						}

						vector<pair<double, const TestedParticleInfo*>> donors;
						donors.reserve(global_particles.size());
						for (size_t gp = 0; gp < global_particles.size(); ++gp) {
							double related_score = -1.0;
							for (size_t up = 0; up < unsolved_paths.size(); ++up) {
								double f = Sim(global_particles[gp].BL, MB[unsolved_paths[up]]);
								if (f > related_score) related_score = f;
							}
							if (related_score >= 0.0) {
								donors.push_back({ related_score, &global_particles[gp] });
							}
						}
						sort(donors.begin(), donors.end(),
							[](const pair<double, const TestedParticleInfo*>& a, const pair<double, const TestedParticleInfo*>& b) {
								return a.first > b.first;
							});

						int donor_cnt = min(ESCAPE_K, (int)donors.size());
						if (donor_cnt > 0) {
							vector<double> fit1, fit2;
							fit1.reserve(popsize);
							fit2.reserve(popsize);
							for (int pi = 0; pi < popsize; ++pi) {
								fit1.push_back(population[pi].fitness);
								fit2.push_back(population2[pi].fitness);
							}
							vector<int> worst1 = argsort(fit1); 
							vector<int> worst2 = argsort(fit2); 

							int replace_cnt = min(ESCAPE_K, popsize);
							for (int rp = 0; rp < replace_cnt; ++rp) {
								int idx1 = worst1[rp];
								const TestedParticleInfo* d1 = donors[rp % donor_cnt].second;
								for (int vd = 0; vd < vardim; ++vd) {
									population[idx1].p[vd] = d1->p_info.p[vd];
									population[idx1].speed[vd] = fabs((float)(rand() % (vmax - vmin + 1)) + vmin);
									perbest[idx1][vd] = population[idx1].p[vd];
								}
								population[idx1].fitness = 0.0;
								perbestfit[idx1] = 0.0;

								int idx2 = worst2[rp];
								const TestedParticleInfo* d2 = donors[(rp + donor_cnt / 2) % donor_cnt].second;
								for (int vd = 0; vd < vardim; ++vd) {
									population2[idx2].p[vd] = d2->p_info.p[vd];
									population2[idx2].speed[vd] = fabs((float)(rand() % (vmax - vmin + 1)) + vmin);
									perbest[idx2][vd] = population2[idx2].p[vd];
								}
								population2[idx2].fitness = 0.0;
								perbestfit[idx2] = 0.0;
							}
							group_escape_stage[group_idx] = 1;
							printf("[INTERVENTION_STAGE1] group=%d, iter=%d, donor_cnt=%d, replace_cnt=%d\n",
								(int)group_idx, iter, donor_cnt, replace_cnt);
							fflush(stdout);
						}
					}
					else {
						
						
						for (int pi = 0; pi < popsize; ++pi) {
							for (int vd = 0; vd < vardim; ++vd) {
								population[pi].p[vd] = fabs((float)(rand() % (xmax - xmin)) + xmin);
								population[pi].speed[vd] = fabs((float)(rand() % (vmax - vmin + 1)) + vmin);
								population2[pi].p[vd] = fabs((float)(rand() % (xmax - xmin)) + xmin);
								population2[pi].speed[vd] = fabs((float)(rand() % (vmax - vmin + 1)) + vmin);
								perbest[pi][vd] = population[pi].p[vd];
							}
							population[pi].fitness = 0.0;
							population2[pi].fitness = 0.0;
							perbestfit[pi] = 0.0;
						}

						globestfit = 0.0;
						gloworstfit = 10000.0;
						const int LOCAL_SAMPLING_MAX = sampsize + selsize * MAXGEN;
						for (int si = sampsize; si < LOCAL_SAMPLING_MAX; ++si) {
							outputing[b][si] = 0.0;
							outputing2[b][si] = 0.0;
							for (int vd = 0; vd < vardim; ++vd) {
								sampling[b][si][vd] = 0.0;
								sampling2[b][si][vd] = 0.0;
							}
						}
						cu_sampling[b] = sampsize;
						update = 1;
						group_escape_stage[group_idx] = 0;
						hist.clear();

						printf("[INTERVENTION_STAGE2] group=%d, iter=%d, sampling_reset=%d\n",
							(int)group_idx, iter, sampsize);
						fflush(stdout);
					}
				}
			}
			group_finished = true;
			for (size_t gm = 0; gm < current_group.size(); ++gm) {
				int path_id = current_group[gm];
				if (path_id >= path_start_index && path_id <= path_end_index && flag[path_id] == 0) {
					group_finished = false;
					break;
				}
			}
			int anchor_done = group_finished ? 1 : 0;
			MPI_Send(&anchor_done, 1, MPI_INT, 0, 33, icomm2);
			if (group_finished)
				break;
			double slope = (maxUp - minUp) / MAXGEN;
			double upFre = MAXGEN - (slope * (MAXGEN - iter));
			if (upFre < minUp) {
				upFre = minUp;
			}
			if ((MAXGEN - iter) % (int)upFre == 0) {
				update = 1;
			}
			else {
				update = 0;
			}
			for (i = 0; i < selsize; i++)
			{
				if (index[abs(2 * popsize - i - 1)] < popsize) {
					population[index[abs(2 * popsize - i - 1)]].fitness = superior_inds[i].fitness;
				}
				else {
					population2[index[abs(2 * popsize - i - 1)] - popsize].fitness = superior_inds[i].fitness;
				}
			}
			k = 0;
			
			const int SAMPLING_MAX = sampsize + selsize * MAXGEN;
			int write_start = cu_sampling[b];
			int write_end = std::min(cu_sampling[b] + selsize, SAMPLING_MAX);
			for (i = write_start; i < write_end; i++)
			{
				outputing[b][i] = superior_inds[k].fitness;
				for (j = 0; j < vardim; j++)
				{
					sampling[b][i][j] = superior_inds[k].p[j];
				}
				if (superior_inds[k].fitness >= 0.6) {
					outputing2[b][i] = superior_inds[k].fitness;
					for (j = 0; j < vardim; j++)
					{
						sampling2[b][i][j] = superior_inds[k].p[j];
					}
				}
				else {
					int sel = (rand() % (selsize - 1 - 0 + 1)) + 0;
					outputing2[b][i] = superior_inds[sel].fitness;
					for (j = 0; j < vardim; j++)
					{
						sampling2[b][i][j] = superior_inds[sel].p[j];
					}
				}
				k++;
			}
			cu_sampling[b] = write_end; 
			loserPop = 0;
			for (i = 0; i < popsize; i++) {
				if (population[i].fitness < loser) {
					loserPop++;
				}
				if (population2[i].fitness < loser) {
					loserPop++;
				}
			}
			if (loserPop >= 2 * popsize * lprop) {
				for (i = 0; i < popsize; i++)
				{
					for (j = 0; j < vardim; j++)
					{
						population[i].p[j] = fabs((float)(rand() % (xmax - xmin)) + xmin);
						population[i].speed[j] = fabs((float)(rand() % (vmax - vmin + 1)) + vmin);
						population2[i].p[j] = fabs((float)(rand() % (xmax - xmin)) + xmin);
						population2[i].speed[j] = fabs((float)(rand() % (vmax - vmin + 1)) + vmin);
					}
				}
				continue;
			}
			else {
				int goLeft = 0, goRight = 0;
				for (i = 1; i < popsize - 1; i++)
				{
					if (population[i].fitness >= Olim) {
						continue;
					}
					for (goLeft = i - 1; goLeft >= 0; goLeft--) {
						if (population[goLeft].fitness >= Olim) {
							personLeft = population[goLeft];
							break;
						}
					}
					if (goLeft < 0) {
						personLeft = population[i - 1];
						for (int j = i - 2; j >= 0; j--) {
							if (population[j].fitness > personLeft.fitness) {
								personLeft = population[j];
							}
						}
					}
					for (goRight = i + 1; goRight < popsize; goRight++) {
						if (population[goRight].fitness >= Olim) {
							personRight = population[goRight];
							break;
						}
					}
					if (goRight == popsize) {
						personRight = population[i + 1];
						for (int j = i + 2; j < popsize; j++) {
							if (population[j].fitness > personRight.fitness) {
								personRight = population[j];
							}
						}
					}
					for (int x = 0; x < vardim; x++) {
						population[i].p[x] = Fgood1 * population[i].p[x] + (1 - Fgood1) * 0.5 * (personLeft.p[x] + personRight.p[x]);
						population[i].speed[x] = Fgood1 * population[i].speed[x] + (1 - Fgood1) * 0.5 * (personLeft.speed[x] + personRight.speed[x]);
					}
				}
				goLeft = 0, goRight = 0;
				for (i = 1; i < popsize - 1; i++)
				{
					if (population2[i].fitness >= Olim) {
						continue;
					}
					for (goLeft = i - 1; goLeft >= 0; goLeft--) {
						if (population2[goLeft].fitness >= Olim) {
							personLeft = population2[goLeft];
							break;
						}
					}
					if (goLeft < 0) {
						personLeft = population2[i - 1];
						for (int j = i - 2; j >= 0; j--) {
							if (population2[j].fitness > personLeft.fitness) {
								personLeft = population2[j];
							}
						}
					}
					for (goRight = i + 1; goRight < popsize; goRight++) {
						if (population2[goRight].fitness >= Olim) {
							personRight = population2[goRight];
							break;
						}
					}
					if (goRight == popsize) {
						personRight = population2[i + 1];
						for (int j = i + 2; j < popsize; j++) {
							if (population2[j].fitness > personRight.fitness) {
								personRight = population2[j];
							}
						}
					}
					for (int x = 0; x < vardim; x++) {
						population2[i].p[x] = Fgood2 * population2[i].p[x] + (1 - Fgood2) * 0.5 * (personLeft.p[x] + personRight.p[x]);
						population2[i].speed[x] = Fgood2 * population2[i].speed[x] + (1 - Fgood2) * 0.5 * (personLeft.speed[x] + personRight.speed[x]);
					}
				}
				estfit2.clear();
				for (i = 0; i < popsize; i++)
				{
					estfit2.push_back(population[i].fitness);
				}
				index2.clear();
				index2 = argsort(estfit2);
				estfit1.clear();
				for (i = 0; i < popsize; i++)
				{
					estfit1.push_back(population2[i].fitness);
				}
				index1.clear();
				index1 = argsort(estfit1);
				for (i = 0; i < NPsQ; i++) {
					for (j = 0; j < vardim; j++)
					{
						population[index2[i]].p[j] = fabs((float)(rand() % (xmax - xmin)) + xmin);
						population[index2[i]].speed[j] = fabs((float)(rand() % (vmax - vmin + 1)) + vmin);
						population2[index1[i]].p[j] = fabs((float)(rand() % (xmax - xmin)) + xmin);
						population2[index1[i]].speed[j] = fabs((float)(rand() % (vmax - vmin + 1)) + vmin);
					}
				}
				for (i = 0; i < popsize; i++)
				{
					if (population[i].fitness > perbestfit[i])
					{
						perbestfit[i] = population[i].fitness;
						for (j = 0; j < vardim; j++)
						{
							perbest[i][j] = population[i].p[j];
						}
					}
					if (perbestfit[i] > globestfit)
					{
						globestfit = perbestfit[i];
						for (j = 0; j < vardim; j++)
						{
							globest[j] = perbest[i][j];
						}
					}
					if (population[i].fitness < gloworstfit)
					{
						gloworstfit = population[i].fitness;
						for (j = 0; j < vardim; j++)
						{
							gloworst[j] = population[i].p[j];
						}
					}
				}
				for (i = 0; i < popsize; i++)
				{
					if (population2[i].fitness > perbestfit[i])
					{
						perbestfit[i] = population2[i].fitness;
						for (j = 0; j < vardim; j++)
						{
							perbest[i][j] = population2[i].p[j];
						}
					}
					if (perbestfit[i] > globestfit)
					{
						globestfit = perbestfit[i];
						for (j = 0; j < vardim; j++)
						{
							globest[j] = perbest[i][j];
						}
					}
					if (population2[i].fitness < gloworstfit)
					{
						gloworstfit = population2[i].fitness;
						for (j = 0; j < vardim; j++)
						{
							gloworst[j] = population2[i].p[j];
						}
					}
				}
				w = wmax - (((wmax - wmin) / MAXGEN) * iter);
				for (i = 0; i < popsize; i++)
				{
					temp = -10000.0;
					for (k = 0; k < popsize; k++)
					{
						temp1 = 0.0;
						temp2 = 0.0;
						temp3 = 0.0;
						temp4 = 0.0;
						for (j = 0; j < vardim; j++)
						{
							temp1 += pow((xmax - xmin), 2);
							temp2 += pow((perbest[k][j] - perbest[i][j]), 2);
						}
						temp3 = sqrt(temp1) / (globestfit - gloworstfit);
						temp4 = sqrt(temp2);
						if ((temp3 * ((perbestfit[k] - perbestfit[i]) / temp4)) > temp)
						{
							temp = temp3 * ((perbestfit[k] - perbestfit[i]));
							tempint = k;
						}
					}
					for (j = 0; j < vardim; j++)
					{
						r1 = (rand() % 1000) / 1000.0;
						r2 = (rand() % 1000) / 1000.0;
						population[i].speed[j] = w * population[i].speed[j] + c1 * r1 * (perbest[i][j] - population[i].p[j]) + c2 * r2 * (perbest[tempint][j] - population[i].p[j]);
						population2[i].speed[j] = w * population2[i].speed[j] + c1 * r1 * (perbest[i][j] - population2[i].p[j]) + c2 * r2 * (perbest[tempint][j] - population2[i].p[j]);
						if (population[i].speed[j] > vmax)
						{
							population[i].speed[j] = vmax;
						}
						if (population[i].speed[j] < vmin)
						{
							population[i].speed[j] = vmin;
						}
						if (population2[i].speed[j] > vmax)
						{
							population2[i].speed[j] = vmax;
						}
						if (population2[i].speed[j] < vmin)
						{
							population2[i].speed[j] = vmin;
						}
					}
					for (j = 0; j < vardim; j++)
					{
						population[i].p[j] = population[i].p[j] + population[i].speed[j];
						population2[i].p[j] = population2[i].p[j] + population2[i].speed[j];
						if (population[i].p[j] > xmax)
						{
							population[i].p[j] = xmax;
						}
						if (population[i].p[j] < xmin)
						{
							population[i].p[j] = xmin;
						}
						if (population2[i].p[j] > xmax)
						{
							population2[i].p[j] = xmax;
						}
						if (population2[i].p[j] < xmin)
						{
							population2[i].p[j] = xmin;
						}
					}
				}
			}
		}
		endtime = MPI_Wtime();
		for (size_t gm = 0; gm < current_group.size(); ++gm) {
			int path_id = current_group[gm];
			if (path_id < path_start_index || path_id > path_end_index) {
				continue;
			}
			if (flag[path_id] == 0) {
				consumedtimes[path_id] = endtime - starttime;
				count[path_id] = iter;
			}
		}
	}
	
	delete[] BL;
	delete[] samples;
	delete[] outputs;
	delete[] samples2;
	delete[] outputs2;
	delete[] sampBL;
	delete[] sampling;
	delete[] outputing;
	delete[] sampling2;
	delete[] outputing2;
	
	int solved_path_count = 0;
	for (b = 0; b < numpath; b++) {
		if (flag[b] == 1) {
			solved_path_count++;
		}
	}
	printf("\n");
	printf("==================================================\n");
	printf("               [EXPERIMENT SUMMARY]\n");
	printf("==================================================\n");
	printf("[Total Time]       : %.4f s\n", MPI_Wtime() - starttime);
	printf("[Successful Paths] : %d/%d\n", solved_path_count, numpath);
	printf("[Executions]       : %d\n", program_executions);
	printf("[Exp Pool Size]    : %zu\n", global_particles.size());
	printf("[Max Gens]         : %d\n", MAXGEN);
	printf("==================================================\n");
	fflush(stdout);
	
	for (b = path_start_index; b <= path_end_index; b++) {
		if (flag[b] == 1) {
			printf("[Path %d] [SUCCESS] - Time: %.4fs, Iterations: %d\n",
				b, consumedtimes[b], count[b]);
			printf("         Solution: ");
			for (j = 0; j < vardim; j++) {
				printf("%d%s", testdata[b][j], (j == vardim - 1) ? "\n" : ", ");
			}
		}
		else {
			const char* failure_reason = (count[b] >= MAXGEN)
				? "Max iterations reached"
				: "Path uncovered";
			printf("[Path %d] [FAILED] - Time: %.4fs, Iterations: %d, Reason: %s\n",
				b, consumedtimes[b], count[b], failure_reason);
		}
		fflush(stdout);
	}
	printf("\n==================================================\n");
	printf("[EXPERIMENT COMPLETE] Executed successfully\n");
	printf("==================================================\n");
	fflush(stdout);
	MPI_Finalize();
	return 0;
}
