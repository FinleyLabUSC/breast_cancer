#ifndef IMMUNE_MODEL_CELL_H
#define IMMUNE_MODEL_CELL_H

#include <array>
#include <vector>
#include <cmath>
#include <random>
#include <string>
#include <iostream>
#include <fstream>
#include "RNG.h"


class Cell{
public:

    // initialization
    static unsigned long cell_counter;  // shared among all cells, keeps track of the number of cells created
    unsigned long unique_cell_ID; // unique, immutable ID per cell, used for tracking

    int runtime_index;                  // a runtime index that is updated every timestep. used in the neighbor list (and other places).


    /*
     * FUNCTIONS
     */

    // initialization
    Cell(std::array<double, 2> loc, std::vector<std::vector<double> > &cellParams, int cellType,
    std::vector<std::string> tCellPhenotypeTrajectory, size_t init_tstamp=0);
    void initialize_Cancer_Cell(std::vector<std::vector<double> > &cellParams, size_t init_tstamp=0);
    void initialize_CD8_Cell(std::vector<std::vector<double> > &cellParams, std::vector<std::string> phenotypeTrajectory, size_t init_tstamp);
    void initialize_CD4_Cell(std::vector<std::vector<double> > &cellParams, size_t init_tstamp=0);
    void initialize_Macrophage_Cell(std::vector<std::vector<double> > &cellParams, size_t init_tstamp=0);
    void initialize_NK_Cell(std::vector<std::vector<double> > &cellParams, size_t init_tstamp=0);
    void initialize_MDSC_Cell(std::vector<std::vector<double> > &cellParams, size_t init_tstamp=0);

    // force functions
    std::array<double, 2> attractiveForce(std::array<double, 2> dx, double otherRadius);
    std::array<double, 2> repulsiveForce(std::array<double, 2> dx, double otherRadius);
    void calculateForces(std::array<double, 2> otherX, double otherRadius, int &otherType);
    void resolveForces(double dt, std::array<double, 2> &tumorCenter, double &necroticRadius, double &necroticForce);
    void resetForces();
    void neighboringCells(std::array<double, 2> otherX, int otherCell_runtime_index);
    void neighboringCancerCells(std::array<double, 2> otherX, int otherCellState);


    // overlap functions
    void calculateOverlap(std::array<double, 2> otherX, double otherRadius);
    void resetOverlap();
    void isCompressed();

    // cell behavior functions
    std::array<double, 3> proliferate(double dt);
    std::array<double, 3> proliferate_v2(double dt);
    void prolifState();
    void prolifState_v2(); // this version is used for the updated Cancer cells that have their own cell cycle clocks
    void inherit(std::vector<double> properties);
    std::vector<double> inheritanceProperties();
    void age(double dt, size_t step_count);

    void migrate(double dt, std::array<double, 2> tumorCenter);
    void migrate_NN(double dt);

    void indirectInteractions(double tstep, size_t step_count);
    void directInteractions(int interactingState, std::array<double, 2> interactingX, std::vector<double> interactionProperties, double tstep);
    std::vector<double> directInteractionProperties(int interactingState, size_t step_count);

    // differentiation
    void differentiate(double dt);

    // cell influences
    void addInfluence(std::array<double, 2> otherX, double otherInfluence, int otherType);
    void addChemotaxis(std::array<double, 2> otherX, double otherInfluence, int otherType);
    void clearInfluence();

    // macrophage
    void macrophage_differentiation(double dt);

    // CD4 specific
    void cd4_differentiation(double dt);

    // CD8 specific
    void cd8_setKillProb(size_t step_count);
    void cd8_pdl1Inhibition(std::array<double, 2> otherX, double otherRadius, double otherpdl1, double dt);

    // cancer specific
    void cancer_dieFromCD8(std::array<double, 2> otherX, double otherRadius, double kp, double dt);
    void cancer_gainPDL1(double dt);
    void cancer_dieFromNK(std::array<double,2> otherX, double otherRadius, double kp, double dt);
    void mutate(int cause, double chemoLevel); // cause = 1 (chemo)    cause = 0 ("natural" mutations)

    // NK specific
    void nk_pdl1Inhibition(std::array<double, 2> otherX, double otherRadius, double otherpdl1, double dt);
    void nk_setKillProb(size_t step_count);

    // MDSC specific
    void mdsc_gainPDL1(double dt);

    // other functions
    double calcDistance(std::array<double, 2> otherX);
    double calcInfDistance(double dist, double xth);
    static double calcNorm(std::array<double, 2> dx);
    static std::array<double, 2> unitVector(std::array<double, 2> v);
    void updateID(int idx);

    double setCellCycleLength(int cellState, int cellType);
    int cellAge;
    void set_cellAge(size_t step_count);

    void updateRunTimeIndex(int index);

    /*
     * PARAMETERS
     */

    // location
    std::array<double, 2> x;

    // physical properties
    double radius;
    bool compressed;
    double currentOverlap;
    std::vector<int> neighbors;
    std::vector<std::array<double, 2> > cancer_neighbors;

    // age, division, and lifespan
    double divProb;
    double divProb_base;
    double deathProb;
    bool canProlif;
    int pTypeStateTransition; // for CD8+ T cells only
    double cellCycleLength;

    double prevDivTime;
    double currDivTime;

    // force properties
    double mu;
    double kc;
    double damping;
    double maxOverlap;
    double rmax;
    std::array<double, 2> currentForces;

    // migration
    double migrationSpeed;
    double migrationBias;

    // cancer properties
    double pdl1Shift;

    // interactions with other cells
    double influenceRadius;
    double pdl1;
    double pdl1WhenExpressed;
    double pdl1_increment;
    std::array<double, 11> influences;
    std::array<double, 8> chemotaxVals;
    double probTh;

    // differentiation
    double kTr;
    double kM1;
    double kM2;
    double plasticity;

    // T cell killing
    double killProb;
    double baseKillProb;
    double infScale;
    std::vector<std::string> t_cell_phenotype_Trajectory;

    std::vector<std::array<double, 2> > location_history;
    void printLocations(int cellNum);

    // chemotherapy
    double chemoDamage;
    double chemoRepair;
    double chemoTolerance;
    double chemoTolRate;
    double deathTime;
    double chemo_Accumulated_Threshold;
    double chemoAccumulated;
    double chemoTime;
    double chemoTimeThresh;
    double chemoUptake;

    double antiPDL1_bindingRate;

    // mutation probability
    double mutationProbability_chemo;
    double mutationProbability_inherent;


    // identification
    int id;
    int type;
    int state;

    //lifespan
    size_t init_time;

    // Cell cycle length
    size_t cellCyclePos;

    size_t birthTime;


private:
    std::mt19937 mt;
};

#endif //IMMUNE_MODEL_CELL_H
