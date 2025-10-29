#include "../inc/RS_Cell.h"
#include <iostream>

// This is a generalized stromal cell that exists only to take up space
Stromal::Stromal(std::array<double, 2> loc, std::vector<std::vector<double>>& cellParams, int cellType, size_t init_tstamp): RS_Cell(loc, cellParams, cellType, init_tstamp)
{
    state = 13;
    // Material properties are inherited from cancer cells
    mu = cellParams[0][0];
    kc = cellParams[1][0];
    damping = cellParams[2][0];
    maxOverlap = cellParams[3][0]*cellParams[4][0];
    radius = cellParams[4][0]/2.0;
    deathProb = 0; // Stromal cells do not die
    migrationSpeed = 0; // Stromal cells do not move
    influenceRadius = cellParams[7][0];
    rmax = 1.5*radius*2;
    init_time = init_tstamp;
    location_history.push_back(x);
}

void Stromal::initialize_cell_from_file(int cell_state, int cell_list_length, double mean_cancer_cell_cycle_length, double std_cancer_cell_cycle_length, RNG& master_rng)
{
    if (cell_state != 13)
    {
        throw std::runtime_error("Stromal cell cannot be initialized with cell_state = " + std::to_string(cell_state) + ".");
    }
    state = cell_state;
    runtime_index = cell_list_length;
    birthTime = 0;
    mother_uniqueID = -1; // A mother_uniqueID of -1 indicates that the cell was part of the initialization of the tumor.
}
