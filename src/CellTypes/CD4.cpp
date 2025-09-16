#include "../inc/RS_Cell.h"

// The child cell class is initialized w/ the parent constructor
// This allows access to specific pre-set variables
CD4::CD4(std::array<double, 2> loc, std::vector<std::vector<double>>& cellParams, int cellType, size_t init_tstamp): RS_Cell(loc, cellParams, cellType, init_tstamp){}

void CD4::initialize_cell(std::vector<std::vector<double>>& cellParams, size_t init_tstamp)
{
    state = 4;
    mu = cellParams[0][1];
    kc = cellParams[1][1];
    damping = cellParams[2][1];
    maxOverlap = cellParams[3][1]*cellParams[4][1];
    radius = cellParams[4][1]/2.0;
    deathProb = cellParams[5][1];
    migrationSpeed = cellParams[6][1];
    kTr = cellParams[7][1];
    influenceRadius = cellParams[8][1];
    migrationBias = cellParams[9][1];
    rmax = 1.5*radius*2;
}

void CD4::initialize_cell_from_file(int cell_state, int cell_list_length, double mean_cancer_cell_cycle_length, double std_cancer_cell_cycle_length, RNG& master_rng)
{
    if (cell_state != 4 && cell_state != 5)
    {
        throw std::runtime_error("CD4 cell cannot be initialized with cell state = " + std::to_string(cell_state) + ".");
    }
    state = cell_state;
    runtime_index = cell_list_length;
    birthTime = master_rng.uniform(0, 1 / deathProb);
    mother_uniqueID = -1; // A mother_uniqueID of -1 indicates that the cell was part of the initialization of the tumor.

    // Only CD4+ Tregs express PDL1; CD4+ Th do not
    if (state == 5)
    {
        pdl1_expression_level = master_rng.uniform(0,max_pdl1_level);
    }
}

