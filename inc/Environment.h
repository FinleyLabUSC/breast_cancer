#ifndef IMMUNE_MODEL_ENVIRONMENT_H
#define IMMUNE_MODEL_ENVIRONMENT_H

#include <vector>
#include <algorithm>
#include <random>
#include "Cell.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <filesystem>
#include <array>


class Environment{

public:
    Environment(std::string folder, std::string set, int base_seed);
    //destructor needed
    void simulate(double tstep, int tx, int met, double binding_rate_pd1_drug);

    void shuffleCells();
    void generateNums();

    // Delete these safely, I don't think they're used.
    std::vector<double> cancer_CellCycle_CDF;
    std::vector<double> CD8_CellCycle_CDF;
    std::vector<double> chemo_treatment_schedule;

    std::vector<double> anti_pd1_treatment_schedule;
    std::vector<double> anti_ctl4_treatment_schedule;

    double effect_of_anti_CTLA4();

    double sampleNormal();

    int model_time;
    void initializeCellsFromFile(std::string filePathway);

private:
    void runCells(double tstep, size_t step_count);
    void neighborInfluenceInteractions(double tstep, size_t step_count);
    void internalCellFunctions(double tstep, size_t step_count);
    void recruitImmuneCells(double tstep, size_t step_count);
    std::array<double, 2> generate_random_location_for_immune_recruitment();
    void tumorSize();
    void necrosis(double tstep);
    double calculateDiffusibles(std::array<double, 2> x);

    void treatment(int tx_flag);

    void anti_pd1(double tstep);
    void anti_pd1_drug(double tstep,double new_dose);
    void anti_ctla4_drug(double tstep,double new_dose);

    void mutateCells();

    void save(double tstep, double tstamp);
    void recordPopulation(double tstamp);
    void record_proliferation(double tstep, int prolifcount);
    void record_drug(double tstep, int tx_type);

    void record_effect(int cellID, double posInfluence, double drug_effect, double ctla4_effect, double scale, double divProb);
    void loadParams();

    void initializeCells();
    void initializeTesting();
    void initializeInVitro();
    void calculateForces(double tstep, size_t step_count);


    void removeDeadCells();
    void updateCell_list();

    void countPops_updateTimeSeries();
    void printStep(double time);

    void neighboringCancerCells(Cell otherCell);
    std::array<std::vector<Cell>,12> createSubLists();

    double nearestNeighborRadius(int state1, int state2);
    
    double dt;

    // cell lists
    std::vector<Cell> cell_list;
    // time courses
    std::vector<int> cancerTS;
    std::vector<int> cd8TS;
    std::vector<int> cd4_th_TS;
    std::vector<int> cd4_treg_TS;
    std::vector<int> m0TS;
    std::vector<int> m1TS;
    std::vector<int> m2TS;
    std::vector<int> nkTS;
    std::vector<int> mdscTS;
    std::vector<double> radiusTS;

    std::vector<double> anti_pd1_TS;
    std::vector<double> anti_ctla4_TS;


    // parameter lists
    std::vector<std::vector<double> > cellParams;
    std::vector<double> recParams;
    std::vector<double> envParams;

    std::string saveDir;
    int steps;
    std::vector<double> immuneCellRecRates;
    std::vector<int> immuneCellRecTypes;
    std::vector<double> immuneCells2rec;
    double recDist;
    double maxRecCytoConc;
    double tumorRadius;
    double necroticGrowth;
    double necroticRadius;
    double necroticForce;
    double necroticLimit;
    std::array<double, 2> tumorCenter;
    double recruitmentDelay;

    // environment params
    double simulationDuration;
    int day;

    // treatment related parameters

    double anti_pd1_decay_rate;
    double anti_ctla4_decay_rate;

    double dose_anti_pd1;
    double dose_anti_ctla4;

    double binding_rate_pd1_drug;

    double mean_cancer_cell_cycle_length = 17;
    double std_cancer_cell_cycle_length = 2;

    RNG rng;

};

#endif //IMMUNE_MODEL_ENVIRONMENT_H
