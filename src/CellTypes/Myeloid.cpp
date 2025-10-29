#include "../inc/RS_Cell.h"
#include <iostream>

// This is a generalized myeloid cell that exists only to take up space
Myeloid::Myeloid(std::array<double, 2> loc, std::vector<std::vector<double>>& cellParams, int cellType, size_t init_tstamp): RS_Cell(loc, cellParams, cellType, init_tstamp)
{
    state = 11;
    // These are the macrophage parameter set minus kM1 and kM2
    mu = cellParams[0][3];
    kc = cellParams[1][3];
    damping = cellParams[2][3];
    maxOverlap = cellParams[3][3]*cellParams[4][3];
    radius = cellParams[4][3]/2.0;
    deathProb = cellParams[5][3]; // Fixed death prob based on macrophages
    migrationSpeed = 0; // Generalized cells do not migrate
    influenceRadius = cellParams[9][3];
    rmax = 1.5*radius*2;
    init_time = init_tstamp;
    location_history.push_back(x);
}

void Myeloid::initialize_cell_from_file(int cell_state, int cell_list_length, double mean_cancer_cell_cycle_length, double std_cancer_cell_cycle_length, RNG& master_rng)
{
    if (cell_state != 11)
    {
        throw std::runtime_error("Myeloid cell cannot be initialized with cell_state = " + std::to_string(cell_state) + ".");
    }
    state = cell_state;
    runtime_index = cell_list_length;
    birthTime = master_rng.uniform(0, 1 / deathProb);
    mother_uniqueID = -1; // A mother_uniqueID of -1 indicates that the cell was part of the initialization of the tumor.
}
