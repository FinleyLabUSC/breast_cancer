//
// Created by Rebecca Bekker on 3/26/25.
//

#include "../inc/RNG.h"

#include <random>
#include <map>
#include "omp.h"
#include <string>
#include <thread>

thread_local std::mt19937 RNG::thread_generator;
thread_local bool RNG::is_initialized = false;

RNG::RNG(int seed) : base_seed(seed) {
    if (!is_initialized) {
        thread_generator.seed(base_seed + std::hash<std::thread::id>{}(std::this_thread::get_id()));
        is_initialized = true;
    }
}

double RNG::uniform(double min, double max) {
    std::uniform_real_distribution<double> dist(min, max);
    return dist(thread_generator);
}

int RNG::uniform_int(int min, int max) {
    std::uniform_int_distribution<int> dist(min, max);
    return dist(thread_generator);
}

double RNG::normal(double mean, double stddev) {
    std::normal_distribution<double> dist(mean, stddev);
    return dist(thread_generator);
}

int RNG::poisson(double mean) {
    std::poisson_distribution<int> dist(mean);
    return dist(thread_generator);
}

std::uniform_real_distribution<double>& RNG::getUniformDistribution(const std::string& name, double min, double max) {
    if (uniformDists.find(name) == uniformDists.end()) {
        uniformDists[name] = std::uniform_real_distribution<double>(min, max);
    }
    return uniformDists[name];
}

std::uniform_int_distribution<int>& RNG::getUniformIntDistribution(const std::string& name, int min, int max){
    if (uniform_Int_Dists.find(name) == uniform_Int_Dists.end()) {
        uniform_Int_Dists[name] = std::uniform_int_distribution<int>(min, max);
    }
    return uniform_Int_Dists[name];
}

std::mt19937& RNG::getGenerator() {
    return thread_generator;
}
