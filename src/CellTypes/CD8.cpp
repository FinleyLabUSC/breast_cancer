#include "../inc/RS_Cell.h"

// The child cell class is initialized w/ the parent constructor
// This allows access to specific pre-set variables
CD8::CD8(std::array<double, 2> loc, std::vector<std::vector<double>>& cellParams, int cellType, size_t init_tstamp): RS_Cell(loc, cellParams, cellType, init_tstamp)
{
    state = 6;
    mu = cellParams[0][2];
    kc = cellParams[1][2];
    damping = cellParams[2][2];
    maxOverlap = cellParams[3][2]*cellParams[4][2];
    radius = cellParams[4][2]/2.0;
    deathProb = cellParams[5][2];
    migrationSpeed = cellParams[6][2];
    next_migrationSpeed = migrationSpeed;
    baseKillProb = cellParams[7][2];
    infScale = cellParams[8][2];
    influenceRadius = cellParams[9][2];
    migrationBias = cellParams[10][2];
    divProb_base = cellParams[11][2];
    divProb = divProb_base;
    deathScale = cellParams[12][2];
    migScale = cellParams[13][2];
    migration_speed_base =migrationSpeed;
    kill_prob_base = baseKillProb;
    killProb = baseKillProb;
    location_history.push_back(x);
    death_prob_base = deathProb;
    rmax = 1.5*radius*2;
    init_time = init_tstamp;
}


void CD8::initialize_cell_from_file(int cell_state, int cell_list_length, double mean_cancer_cell_cycle_length, double std_cancer_cell_cycle_length, RNG& master_rng)
{
    if (cell_state != 6)
    {
        throw std::runtime_error("CD8 cell cannot be initialized with cell state = " + std::to_string(cell_state) + ".");
    }
    state = cell_state;
    runtime_index = cell_list_length;
    birthTime = master_rng.uniform(0, 1 / deathProb);
    mother_uniqueID = -1; // A mother_uniqueID of -1 indicates that the cell was part of the initialization of the tumor.

    // CD8-specific properties -- we assume cells are exhausted in mIHC imaging
    pd1_expression_level = master_rng.uniform(0,max_pd1_level);
    deathProb = master_rng.uniform(4 * death_prob_base, 20 * death_prob_base); // High death rate
    migrationSpeed = master_rng.uniform(0, 0.25 * migration_speed_base); // Low migration speed
    divProb = master_rng.uniform(0, 0.25 * divProb_base); // Low division probability
    killProb = master_rng.uniform(0, 0.25 * kill_prob_base); // Low cyctotoxic effect
}
