// variant_a.h
#pragma once
#include <vector>
#include <algorithm>
#include <random>
#include <ctime>
#include <iostream>
#include "variable.h"


void variantA_initialize(individual* population,
    individual* population2,
    const double samples[][vardim],
    const double outputs[][sampsize],
    int path_index);