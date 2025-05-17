//
// Created by Rebecca Bekker on 3/26/25.
//

#ifndef RNG_H
#define RNG_H

#include <random>
#include <map>
#include <string>

class RNG {
public:
    RNG(int seed);
    double uniform(double min, double max);
    int uniform_int(int min, int max);
    double normal(double mean, double stddev);
    int poisson(double mean);

    // Optional: Methods for creating distributions on demand
    std::uniform_real_distribution<double>& getUniformDistribution(const std::string& name, double min, double max);
    std::uniform_int_distribution<int>& getUniformIntDistribution(const std::string& name, int min, int max);
    std::normal_distribution<double>& getNormalDistribution(const std::string& name, double mean, double stddev);
    std::poisson_distribution<int>& getPoissonDistribution(const std::string& name, double mean);

    std::mt19937& getGenerator();

private:
    int base_seed;
    static thread_local std::mt19937 thread_generator;
    static thread_local bool is_initialized;

    //std::mt19937 generator;

    // Optional: Store distributions for reuse (if needed)
    std::map<std::string, std::uniform_real_distribution<double>> uniformDists;
    std::map<std::string, std::uniform_int_distribution<int>> uniform_Int_Dists;
    std::map<std::string, std::normal_distribution<double>> normalDists;
    std::map<std::string, std::poisson_distribution<int>> poissonDists;
};



#endif //RNG_H
