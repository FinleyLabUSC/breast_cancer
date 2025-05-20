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
    Environment(std::string folder, std::string set, std::string tCellTrajectoryPath);
    //destructor needed
    void simulate(double tstep, int tx, int met, int numDoses);
    void truncatedGaussian(double mean, double stddev, double min, double max);

    // Visualization on / off
    bool visualize;

    // Delete these safely, I don't think they're used.
    std::vector<double> cancer_CellCycle_CDF;
    std::vector<double> CD8_CellCycle_CDF;
    std::vector<double> chemo_treatment_schedule;

    std::vector<double> ici_treatment_schedule;

    double sampleNormal();

    int model_time;
    void initializeCellsFromFile(std::string filePathway);

private:
    void runCells(double tstep, size_t step_count);
    void neighborInfluenceInteractions(double tstep, size_t step_count);
    void internalCellFunctions(double tstep, size_t step_count);
    void recruitImmuneCells(double tstep, size_t step_count);
    void recruitImmuneCells_random(double tstep, size_t step_count);
    std::array<double, 2> recruitmentLocation();
    std::array<double, 2> recruitmentLocation_random();
    void tumorSize();
    void necrosis(double tstep);
    double calculateDiffusibles(std::array<double, 2> x);

    void treatment();
    void chemotherapy(double tstep);
    void chemotherapy_drug(double tstep, double new_dose);

    void immune_checkpoint_inhibitor(double tstep);
    void immune_checkpoint_inhibitor_drug(double tstep,double new_dose);


    void mutateCells(int cause);

    void populateTrajectories(std::string tCellTrajectoryPath);
    std::vector<std::string> getTcellTrajectory(); // This function will check whether the tCellPhenotypeTrajectory's have been read in, read the files if not, and return a randomly selected trajectory

    void save(double tstep, double tstamp);
    void recordPopulation(double tstamp);
    void loadParams();

    void initializeCells();
    void initializeTesting();
    void initializeInVitro();
    void calculateForces(double tstep);


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
    std::vector<double> chemoTS;
    std::vector<double> ICI_TS;
    /*
    t cell trajectory matrix: we can either represent as a vector of chars 
    where a char maps to a phenotypic state or an int where the int maps to 
    a phenotypic state
    */
    std::vector<std::string> tCellPhenotypeTrajectory_1;

    std::vector<std::vector<std::string> > tCellPhenotypeTrajectory;

    std::string tCellTrajectoryPathway;

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

    double chemotherapy_decay_rate;
    double ICI_decay_rate;
    double dose_chemo;
    double dose_ICI;
    int number_of_chemo_doses;

    std::mt19937 mt;

};

#endif //IMMUNE_MODEL_ENVIRONMENT_H
