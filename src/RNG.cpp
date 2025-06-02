//
// Created by Rebecca Bekker on 3/26/25.
//

#include "../inc/RNG.h"

#include <random>
#include <map>
#include "omp.h"
#include <string>
#include <thread>


RNG::RNG(int seed) : base_seed(seed) {
    instance_generator.seed(seed);
}

std::mt19937 &RNG::getGenerator() {
    return instance_generator;
}

unsigned int RNG::get_context_seed(size_t timestep, unsigned long uniqueCellID, int processID) {
    unsigned int h = base_seed;
    h ^= timestep + 0x9e3779b9 + (h << 6) + (h >> 2); // Fibonacci hash constant
    h ^= uniqueCellID + 0x9e3779b3 + (h << 6) + (h >> 2); // Another constant
    h ^= processID + 0x9e3779b1 + (h << 6) + (h >> 2); // Another constant
    return h;
}

double RNG::uniform(double min, double max) {
    std::uniform_real_distribution<double> dist(min, max);
    return dist(instance_generator);
}

int RNG::uniform_int(int min, int max) {
    std::uniform_int_distribution<int> dist(min, max);
    return dist(instance_generator);
}

double RNG::normal(double mean, double stddev) {
    std::normal_distribution<double> dist(mean, stddev);
    return dist(instance_generator);
}

int RNG::poisson(double mean) {
    std::poisson_distribution<int> dist(mean);
    return dist(instance_generator);
}


double RNG::uniform(double min, double max, std::mt19937& gen) {
    std::uniform_real_distribution<double> dist(min, max);
    return dist(gen);
}

int RNG::uniform_int(int min, int max, std::mt19937& gen) {
    std::uniform_int_distribution<int> dist(min, max);
    return dist(gen);
}

double RNG::normal(double mean, double stddev, std::mt19937& gen) {
    std::normal_distribution<double> dist(mean, stddev);
    return dist(gen);
}

int RNG::poisson(double mean, std::mt19937& gen) {
    std::poisson_distribution<int> dist(mean);
    return dist(gen);
}
