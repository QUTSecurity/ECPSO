// PopulationInit.cpp
#pragma once
#include <numeric>
#include <random>
#include "variable.h"

// Initialize populations with fitness-guided sampling and random injection.
void FitnessGuidedInitialize(individual* population,
    individual* population2,
    const double samples[][vardim],
    const double outputs[][sampsize],
    int pathIndex) {

    const double alpha = 1;

    double base_delta = (xmax - xmin) * 0.5;
    double scale = sqrt((double)vardim / sampsize);
    double delta = base_delta * scale / vardim;

    //const double delta = (xmax - xmin) * (0.5 / (double)vardim);

    const double rand_ratio = 0.20;
    const int rand_num = static_cast<int>(popsize * rand_ratio);
    const int guide_num = popsize - rand_num;

    std::vector<double> fitnesses(sampsize);
    for (int i = 0; i < sampsize; ++i) {
        fitnesses[i] = outputs[pathIndex][i];
    }

    std::vector<int> sorted_idx(sampsize);
    std::iota(sorted_idx.begin(), sorted_idx.end(), 0);
    std::sort(sorted_idx.begin(), sorted_idx.end(), [&fitnesses](int a, int b) {
        return fitnesses[a] > fitnesses[b];
        });
    std::vector<double> weights(sampsize);
    double total = 0.0;
    for (int i = 0; i < sampsize; ++i) {
        weights[i] = (fitnesses[i] > 0 ? fitnesses[i] : 1e-6);
        total += weights[i];
    }
    std::vector<double> cumulative(sampsize);
    cumulative[0] = weights[0] / total;
    for (int i = 1; i < sampsize; ++i) {
        cumulative[i] = cumulative[i - 1] + weights[i] / total;
    }
    std::mt19937 rng(static_cast<unsigned int>(time(nullptr)));
    std::uniform_real_distribution<double> perturb(-delta, delta);
    std::uniform_real_distribution<double> rand01(0.0, 1.0);
    std::uniform_real_distribution<double> speed_dist(vmin, vmax);
    std::uniform_real_distribution<double> val_dist(xmin, xmax);

    int pid = 0;
    // Generate guided particles around sampled regions.
    for (int i = 0; i < guide_num && pid < popsize; ++i) {
        double r = rand01(rng);
        int idx = std::lower_bound(cumulative.begin(), cumulative.end(), r) - cumulative.begin();
        if (idx >= sampsize) idx = sampsize - 1;

        for (int d = 0; d < vardim; ++d) {
            double val = samples[idx][d] + perturb(rng);
            val = std::min(std::max(val, (double)xmin), (double)xmax);
            population[pid].p[d] = val;
            population2[pid].p[d] = val;
            population[pid].speed[d] = speed_dist(rng);
            population2[pid].speed[d] = speed_dist(rng);
        }
        population[pid].fitness = 0.0;
        population2[pid].fitness = 0.0;
        ++pid;
    }
    // Fill the remaining particles with random initialization.
    for (; pid < popsize; ++pid) {
        for (int d = 0; d < vardim; ++d) {
            double val = val_dist(rng);
            population[pid].p[d] = val;
            population2[pid].p[d] = val;
            population[pid].speed[d] = speed_dist(rng);
            population2[pid].speed[d] = speed_dist(rng);
        }
        population[pid].fitness = 0.0;
        population2[pid].fitness = 0.0;
    }
}
