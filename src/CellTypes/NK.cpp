#include "../inc/RS_Cell.h"

// The child cell class is initialized w/ the parent constructor
// This allows access to specific pre-set variables
NK::NK(std::array<double, 2> loc, std::vector<std::vector<double>>& cellParams, int cellType, size_t init_tstamp): RS_Cell(loc, cellParams, cellType, init_tstamp)
{
    state = 8;
    mu = cellParams[0][4];
    kc = cellParams[1][4];
    damping = cellParams[2][4];
    maxOverlap = cellParams[3][4]*cellParams[4][4];
    radius = cellParams[4][4]/2.0;
    deathProb = cellParams[5][4];
    migrationSpeed = cellParams[6][4];
    baseKillProb = cellParams[7][4];
    killProb = baseKillProb;
    infScale = cellParams[8][4];
    influenceRadius = cellParams[9][4];
    migrationBias = cellParams[10][4];
    divProb_base = 0;
    divProb_base = 0;
    deathScale = cellParams[12][4];
    migScale = cellParams[13][4];
    rmax = 1.5*radius*2;
}

void NK::initialize_cell_from_file(int cell_state, int cell_list_length, double mean_cancer_cell_cycle_length, double std_cancer_cell_cycle_length, RNG& master_rng)
{
    if (cell_state != 8)
    {
        throw std::runtime_error("NK cell cannot be initialized with cell_state =" + std::to_string(cell_state) + ".");
    }
    state = cell_state;
    runtime_index = cell_list_length;
    birthTime = master_rng.uniform(0, 1 / deathProb);
    mother_uniqueID = -1; // A mother_uniqueID of -1 indicates that the cell was part of the initialization of the tumor.

    // NK Cells express PD1 & are cytotoxic
    pd1_expression_level = master_rng.uniform(0,max_pd1_level);
    deathProb = master_rng.uniform(4 * death_prob_base, 20 * death_prob_base); // High death rate
    migrationSpeed = master_rng.uniform(0, 0.25 * migration_speed_base);
    killProb = master_rng.uniform(0, 0.25 * kill_prob_base); // Low cytotoxic effect
}

void NK::indirectInteractions(double tstep, size_t step_count, RNG& master_rng, std::mt19937& temporary_rng, double anti_pd1_concentration, double binding_rate_pd1_drug)
{
    // TODO: Update PD1 expression function & update_indirectProperties to match declaration
    express_PD1();
    update_indirectProperties();
}

void NK::directInteractions(int interactingState, std::array<double, 2> interactingX, std::vector<double> interactionProperties, double tstep, RNG& master_gen, std::mt19937& temporary_rng)
{
    if (interactingState == 2 || interactingState == 3 || interactingState == 5 || interactingState == 10)
    {
        // TODO: Update PDL1 inhibition function when implemented
        pdl1_inhibition();
    }
}
