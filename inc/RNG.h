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
    explicit RNG(int seed);

    double uniform(double min, double max);
    int uniform_int(int min, int max);
    double normal(double mean, double stddev);
    int poisson(double mean);

    double uniform(double min, double max, std::mt19937& gen);
    int uniform_int(int min, int max, std::mt19937& gen);
    double normal(double mean, double stddev, std::mt19937& gen);
    int poisson(double mean, std::mt19937& gen);

    unsigned int get_context_seed(size_t timestep, unsigned long uniqueCellIS, int processID);
    std::mt19937& getGenerator();
private:
    int base_seed;
    std::mt19937 instance_generator;
};



#endif //RNG_H
