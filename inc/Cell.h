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
    Cell(std::array<double, 2> loc, std::vector<std::vector<double> > &cellParams, int cellType, size_t init_tstamp=0);
    void initialize_cell_from_file(int state, int run_time_index, double mean_cancer_cell_cycle_length, double std_cancer_cell_cycle_length, RNG& master_rng);
    void initialize_Cancer_Cell(std::vector<std::vector<double> > &cellParams, size_t init_tstamp=0);
    void initialize_CD8_Cell(std::vector<std::vector<double> > &cellParams, size_t init_tstamp);
    void initialize_CD4_Cell(std::vector<std::vector<double> > &cellParams, size_t init_tstamp=0);
    void initialize_Macrophage_Cell(std::vector<std::vector<double> > &cellParams, size_t init_tstamp=0);
    void initialize_NK_Cell(std::vector<std::vector<double> > &cellParams, size_t init_tstamp=0);
    void initialize_MDSC_Cell(std::vector<std::vector<double> > &cellParams, size_t init_tstamp=0);

    // force functions
    std::array<double, 2> attractiveForce(std::array<double, 2> dx, double otherRadius);
    std::array<double, 2> repulsiveForce(std::array<double, 2> dx, double otherRadius);
    void calculateForces(std::array<double, 2> otherX, double otherRadius, int &otherType);
    void resolveForces(double dt, RNG& master_rng, std::mt19937& temporary_rng);
    void resetForces(RNG& master_rng, std::mt19937& local_gen);

    void determine_neighboringCells(std::array<double,2> otherX, int otherCell_runtime_index, int otherCell_state);

    // overlap functions
    void calculateOverlap(std::array<double, 2> otherX, double otherRadius);
    void resetOverlap();
    void isCompressed();

    // cell behavior functions
    std::array<double, 3> proliferate(double dt, RNG& master_rng);
    void proliferationState(double anti_ctla4_concentration); // this version is used for the updated Cancer cells that have their own cell cycle clocks
    void inherit(std::vector<double> properties);
    std::vector<double> inheritanceProperties();
    void age(double dt, size_t step_count, RNG& master_rng);

    void migrate(double dt, std::array<double, 2> tumorCenter);
    void migrate_NN(double dt,RNG& master_rng, std::mt19937& temporary_rng);

    void indirectInteractions(double tstep, size_t step_count, RNG& master_rng, std::mt19937& temporary_rng, double anti_pd1_concentration, double binding_rate_pd1_drug);
    void directInteractions(int interactingState, std::array<double, 2> interactingX, std::vector<double> interactionProperties, double tstep, RNG& master_gen, std::mt19937& temporary_rng);
    std::vector<double> directInteractionProperties(int interactingState, size_t step_count);

    // differentiation
    void differentiate(double dt, RNG& master_rng, std::mt19937& temporary_rng);

    // cell influences
    void addInfluence(std::array<double, 2> otherX, double otherInfluence, int otherType);
    void addChemotaxis(std::array<double, 2> otherX, double otherInfluence, int otherType);
    void clearInfluence();

    // macrophage
    void macrophage_differentiation(double dt, RNG& master_rng, std::mt19937& temporary_rng);
    void macrophage_pdl1_expression_level(double dt);

    // CD4 specific
    void cd4_differentiation(double dt, RNG& master_rng, std::mt19937& localgen);
    void cd4_pdl1_expression_level(double dt);

    // CD8 specific
    void cd8_setKillProb(size_t step_count);
    void cd8_pdl1Inhibition(std::array<double, 2> otherX, double otherRadius, double otherpdl1, double dt, RNG& master_rng, std::mt19937& local_gen);
    double cd8_setProliferationScale(double anti_ctla4_concentration);
    void cd8_pd1_expression_level(double dt, double anti_pd1_concentration, double binding_rate_pd1_drug);
    void cd8_update_properties_indirect();

    // cancer specific
    void cancer_dieFromCD8(std::array<double, 2> otherX, double otherRadius, double immune_cell_kill_prob, double dt, RNG& master_rng, std::mt19937& localgen);
    void cancer_gainPDL1(double dt);
    void cancer_dieFromNK(std::array<double,2> otherX, double otherRadius, double immune_cell_kill_prob, double dt, RNG& master_rng, std::mt19937& local_gen);
    void mutate(RNG& master_rng); // cause = 0 ("natural" mutations)

    // NK specific
    void nk_pdl1Inhibition(std::array<double, 2> otherX, double otherRadius, double otherpdl1, double dt,RNG& master_rng, std::mt19937& local_gen);
    void nk_update_properties_indirect(size_t step_count);
    void nk_pd1_expression_level(double dt, double anti_pd1_concentration, double binding_rate_pd1_drug);

    double sensitivity_to_antiPD1(double anti_pd1_concentration, double binding_rate_pd1_drug);
    double sensitivity_to_antiCTLA4();
    // MDSC specific
    void mdsc_gainPDL1(double dt);

    // other functions
    double calcDistance(std::array<double, 2> otherX);
    double calcInfDistance(double dist, double xth);
    static double calcNorm(std::array<double, 2> dx);
    static std::array<double, 2> unitVector(std::array<double, 2> v);
    void updateID(int idx);
    double Hill_function(double concentration, double EC50, double n);

    void set_cell_cycle_length(double cell_cycle_length);

    void set_cellAge(size_t step_count);

    void updateRunTimeIndex(int index);
    void resetImmuneSynapse();

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
    std::vector<std::array<double, 2>> cancer_neighbors;

    int cellAge;

    // age, division, and lifespan
    double divProb;
    double divProb_base;

    double deathProb;
    double death_prob_base;


    double kill_prob_base;

    bool canProlif;

    double cellCycleLength;

    double prevDivTime;
    double currDivTime;

    // force properties
    double mu;
    double kc;
    double damping;
    double maxOverlap;
    double rmax;

    double maxRepulsiveForce = 30;
    std::array<double, 2> currentForces;

    // migration
    double migrationSpeed;
    double migrationBias;
    double migration_speed_base;

    bool immuneSynapseFormed = false; // by default all cells have no immune synapses formed.
    // interactions with other cells
    double influenceRadius;

    // cancer properties
    double pdl1_expression_level;
    double threshold_for_pdl1_induction;
    double pdl1_induction_rate;
    double pdl1_decay;
    double max_pdl1_level;

    double pd1_induction_rate;
    double pd1_drug_bound;
    double pd1_available;
    double pd1_expression_level;
    double threshold_for_pd1_induction;
    double pd1_decay_rate;
    double max_pd1_level;

    double inhibitory_effect_of_binding_PD1_PDL1;// This specifies the level of reduction in kill prob and migration prob from PD1 - PDL1 binding.
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
    double migScale;
    double deathScale;

    // Used to ensure that no data races occur during the parallel regions in runCells
    double next_killProb;
    int next_state;
    double next_migrationSpeed;
    double next_death_prob;

    // anti CTTLA-4 parameters
    double anti_CTLA4_IC50 = 5; // 5 nM ref: https://bpsbioscience.com/ctla4-neutralizing-antibody-71212
    double anti_CTLA4_hill_coeff = 1; // steepness of sigmoid

    // anti PD1 parameters
    double anti_pdl1_hill_coeff = 1; // steepness of sigmoid

    std::vector<std::array<double, 2> > location_history;
    void printLocations(std::string saveDir);

    // mutation probability
    double mutationProbability_inherent;
    double mutation_prob_PDL1;

    // identification
    int id;
    int type;
    int state;
    int mother_uniqueID;

    int death_type = -1;

    //lifespan
    size_t init_time;

    // Cell cycle length
    size_t cellCyclePos;

    size_t birthTime;
};

#endif //IMMUNE_MODEL_CELL_H
