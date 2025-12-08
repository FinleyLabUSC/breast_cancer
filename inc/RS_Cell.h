#ifndef BREAST_CANCER_RS_CELL_H
#define BREAST_CANCER_RS_CELL_H

#include <array>
#include <vector>
#include <cmath>
#include <random>
#include <string>
#include "RNG.h"
#include <iostream>
#include <fstream>

class RS_Cell
{
    public:
    // virtual destructor
    virtual ~RS_Cell() = default;

    // identification
    static unsigned long cell_counter; // counter shared amongst all cells to track # of cells created
    unsigned long unique_cell_ID; // unique immutable ID for each cell
    int runtime_index; // mutable index in cell list updated every step

    // initialization
    RS_Cell(std::array<double, 2> loc, std::vector<std::vector<double>> &cellParams, int cellType, size_t init_tstamp=0);
    virtual void initialize_cell_from_file(int state, int cell_list_length, double mean_cancer_cell_cycle_length, double std_cancer_cell_cycle_length, RNG& master_rng);

    // force functions
    std::array<double, 2> attractiveForce(std::array<double, 2> dx, double otherRadius);
    std::array<double, 2> repulsiveForce(std::array<double, 2> dx, double otherRadius);
    std::array<double, 2> synapse_springForce(std::array<double, 2> dx, double otherRadius, int &otherType);
    void calculateForces(std::array<double, 2> otherX, double otherRadius, int &otherType);
    void resolveForces(double dt, RNG& master_rng, std::mt19937& temporary_rng);
    void resetForces(RNG& master_rng, std::mt19937& local_gen);

    // neighboring cells
    std::array<double, 2> determine_grid();
    void determine_neighboringCells(std::array<double,2> otherX, int otherCell_runtime_index, int otherCell_state);
    void determine_immuneSynapses(std::array<double, 2> otherX, int otherRadius, int& otherType, unsigned long otherUID);

    // overlap functions
    void calculateOverlap(std::array<double, 2> otherX, double otherRadius);
    void resetOverlap();
    void isCompressed();

    // proliferation and aging
    virtual std::array<double, 3> proliferate(double dt, RNG& master_rng);
    virtual void mutate(RNG& master_rng);
    virtual void proliferationState(double anti_ctla4_concentration);
    virtual void inherit(std::vector<double> properties);
    virtual std::vector<double> inheritanceProperties();
    void age(double dt, size_t step_count, RNG& master_rng);
    std::array<double, 3> cycle_proliferate(double dt, RNG& master_rng); // Proliferate on a clock
    std::array<double, 3> prob_proliferate(double dt, RNG& master_rng); // Proliferate with a probability

    // movement
    virtual void migrate_NN(double dt,RNG& master_rng, std::mt19937& temporary_rng);

    // cell-cell interactions
    virtual void indirectInteractions(double tstep, size_t step_count, RNG& master_rng, std::mt19937& temporary_rng, double anti_pd1_concentration, double binding_rate_pd1_drug);
    virtual void directInteractions(int interactingState, std::array<double, 2> interactingX, std::vector<double> interactionProperties, double tstep, RNG& master_gen, std::mt19937& temporary_rng);
    virtual std::vector<double> directInteractionProperties(int interactingState, size_t step_count);
    virtual void express_PDL1(double dt);
    virtual void express_PD1(double dt, double anti_pd1_concentration, double binding_rate_pd1_drug);
    virtual void update_indirectProperties(size_t step_count);
    virtual void contact_die(int killer_state, std::array<double, 2> otherX, double otherRadius, double kill_prob, double dt, RNG& master_rng, std::mt19937& temporary_rng);
    virtual void pdl1_inhibition(std::array<double, 2> otherX, double otherRadius, double otherPDL1, double dt, RNG& master_rng, std::mt19937&
                                 temporary_rng);
    void resetImmuneSynapse();

    // differentiation
    virtual void differentiate(double dt, RNG& master_rng, std::mt19937& temporary_rng);

    // cell influences
    void addInfluence(std::array<double, 2> otherX, double otherInfluence, int otherType);
    void clearInfluence();

    // math functions
    double calcDistance(std::array<double, 2> otherX);
    double calcInfDistance(double dist, double xth);
    static double calcNorm(std::array<double, 2> dx);
    static std::array<double, 2> unitVector(std::array<double, 2> v);
    double Hill_function(double concentration, double EC50, double n);

    // cell cycle functions
    void set_cell_cycle_length(double cell_cycle_length);
    double sensitivity_to_antiPD1(double anti_pd1_concentration, double binding_rate_pd1_drug);
    double sensitivity_to_antiCTLA4();
    void set_cellAge(size_t step_count);

    // update functions
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
    std::vector<int> neighbors; // Vector to store CellIDs of neighbors
    std::vector<unsigned long> synapses; // Vector to store CellIDs of synapses
    std::vector<int> synapse_durations; // Vector to store duration of synapses
    std::vector<int> synapse_types; // Vector to store TYPE of synapse (all Type I right now)
    std::vector<std::array<double, 2>> cancer_neighbors;

    // age, division, and lifespan
    int cellAge;
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

/* CHILD CELL CLASSES:
 * When adding cells to the model, declare them here!
 * Because cells have different properties, each MUST have their own constructor.
 * Currently, cell types are final, meaning they should not have children.
 * To avoid changing types mid-simulation, differentiation should be handled within the class, using the "state" variable.
 */

class Cancer final : public RS_Cell
{
    public:
    Cancer(std::array<double, 2> loc, std::vector<std::vector<double>> &cellParams, int cellType, size_t init_tstamp=0);
    void initialize_cell_from_file(int cell_state, int cell_list_length, double mean_cancer_cell_cycle_length, double std_cancer_cell_cycle_length, RNG& master_rng) override;
    std::array<double, 3> proliferate(double dt, RNG& master_rng) override;
    void migrate_NN(double dt, RNG& master_rng, std::mt19937& temporary_rng) override;
    void indirectInteractions(double tstep, size_t step_count, RNG& master_rng, std::mt19937& temporary_rng, double anti_pd1_concentration, double binding_rate_pd1_drug) override;
    void directInteractions(int interactingState, std::array<double, 2> interactingX, std::vector<double> interactionProperties, double tstep, RNG& master_gen, std::mt19937& temporary_rng) override;
    void contact_die(int killer_state, std::array<double, 2> otherX, double otherRadius, double kill_prob, double dt, RNG& master_rng, std::
                     mt19937& temporary_rng) override;
    std::vector<double> directInteractionProperties(int interactingState, size_t step_count) override;
    void mutate(RNG& master_rng) override;
    void inherit(std::vector<double> properties) override;
    std::vector<double> inheritanceProperties() override;
    void proliferationState(double anti_ctla4_concentration) override;
};

class CD4 final : public RS_Cell
{
    public:
    CD4(std::array<double, 2> loc, std::vector<std::vector<double>> &cellParams, int cellType, size_t init_tstamp=0);
    void initialize_cell_from_file(int cell_state, int cell_list_length, double mean_cancer_cell_cycle_length, double std_cancer_cell_cycle_length, RNG& master_rng) override;
    void indirectInteractions(double tstep, size_t step_count, RNG& master_rng, std::mt19937& temporary_rng, double anti_pd1_concentration, double binding_rate_pd1_drug) override;
    std::vector<double> directInteractionProperties(int interactingState, size_t step_count) override;
    void differentiate(double dt, RNG& master_rng, std::mt19937& temporary_rng) override;
};

class CD8 final : public RS_Cell
{
    public:
    CD8(std::array<double, 2> loc, std::vector<std::vector<double>> &cellParams, int cellType, size_t init_tstamp=0);
    void initialize_cell_from_file(int cell_state, int cell_list_length, double mean_cancer_cell_cycle_length, double std_cancer_cell_cycle_length, RNG& master_rng) override;
    std::array<double, 3> proliferate(double dt, RNG& master_rng) override;
    void indirectInteractions(double tstep, size_t step_count, RNG& master_rng, std::mt19937& temporary_rng, double anti_pd1_concentration, double binding_rate_pd1_drug) override;
    void directInteractions(int interactingState, std::array<double, 2> interactingX, std::vector<double> interactionProperties, double tstep, RNG& master_gen, std::mt19937& temporary_rng) override;
    std::vector<double> directInteractionProperties(int interactingState, size_t step_count) override;
    void update_indirectProperties(size_t step_count) override;
    void inherit(std::vector<double> properties) override;
    std::vector<double> inheritanceProperties() override;
    void proliferationState(double anti_ctla4_concentration) override;
};

class Macrophage final : public RS_Cell
{
    public:
    Macrophage(std::array<double, 2> loc, std::vector<std::vector<double>> &cellParams, int cellType, size_t init_tstamp=0);
    void initialize_cell_from_file(int cell_state, int cell_list_length, double mean_cancer_cell_cycle_length, double std_cancer_cell_cycle_length, RNG& master_rng) override;
    void indirectInteractions(double tstep, size_t step_count, RNG& master_rng, std::mt19937& temporary_rng, double anti_pd1_concentration, double binding_rate_pd1_drug) override;
    std::vector<double> directInteractionProperties(int interactingState, size_t step_count) override;
    void differentiate(double dt, RNG& master_rng, std::mt19937& temporary_rng) override;
};

class MDSC final : public RS_Cell
{
    public:
    MDSC(std::array<double, 2> loc, std::vector<std::vector<double>> &cellParams, int cellType, size_t init_tstamp=0);
    void initialize_cell_from_file(int cell_state, int cell_list_length, double mean_cancer_cell_cycle_length, double std_cancer_cell_cycle_length, RNG& master_rng) override;
    void indirectInteractions(double tstep, size_t step_count, RNG& master_rng, std::mt19937& temporary_rng, double anti_pd1_concentration, double binding_rate_pd1_drug) override;
    std::vector<double> directInteractionProperties(int interactingState, size_t step_count) override;
};

class NK final : public RS_Cell
{
    public:
    NK(std::array<double, 2> loc, std::vector<std::vector<double>> &cellParams, int cellType, size_t init_tstamp=0);
    void initialize_cell_from_file(int cell_state, int cell_list_length, double mean_cancer_cell_cycle_length, double std_cancer_cell_cycle_length, RNG& master_rng) override;
    void indirectInteractions(double tstep, size_t step_count, RNG& master_rng, std::mt19937& temporary_rng, double anti_pd1_concentration, double binding_rate_pd1_drug) override;
    void directInteractions(int interactingState, std::array<double, 2> interactingX, std::vector<double> interactionProperties, double tstep, RNG& master_gen, std::mt19937& temporary_rng) override;
    std::vector<double> directInteractionProperties(int interactingState, size_t step_count) override;
    void update_indirectProperties(size_t step_count) override;
};

class Lymphoid final : public RS_Cell
{
    public:
    Lymphoid(std::array<double, 2> loc, std::vector<std::vector<double>> &cellParams, int cellType, size_t init_tstamp=0);
    void initialize_cell_from_file(int cell_state, int cell_list_length, double mean_cancer_cell_cycle_length, double std_cancer_cell_cycle_length, RNG& master_rng) override;
    void inherit(std::vector<double> properties) override;
    std::vector<double> inheritanceProperties() override;
    std::array<double, 3> proliferate(double dt, RNG& master_rng) override;
};

class Myeloid final : public RS_Cell
{
    public:
    Myeloid(std::array<double, 2> loc, std::vector<std::vector<double>> &cellParams, int cellType, size_t init_tstamp=0);
    void initialize_cell_from_file(int cell_state, int cell_list_length, double mean_cancer_cell_cycle_length, double std_cancer_cell_cycle_length, RNG& master_rng) override;
};

class Stromal final : public RS_Cell
{
    public:
    Stromal(std::array<double, 2> loc, std::vector<std::vector<double>> &cellParams, int cellType, size_t init_tstamp=0);
    void initialize_cell_from_file(int cell_state, int cell_list_length, double mean_cancer_cell_cycle_length, double std_cancer_cell_cycle_length, RNG& master_rng) override;
};

#endif //BREAST_CANCER_RS_CELL_H