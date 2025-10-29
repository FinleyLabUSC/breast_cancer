#include "../inc/RS_Cell.h"
#include <iostream>

// This is a generalized lymphoid cell that exists only to take up space
Lymphoid::Lymphoid(std::array<double, 2> loc, std::vector<std::vector<double>>& cellParams, int cellType, size_t init_tstamp): RS_Cell(loc, cellParams, cellType, init_tstamp)
{
    state = 12;
    // This is the CD8 parameter set
    mu = cellParams[0][2];
    kc = cellParams[1][2];
    damping = cellParams[2][2];
    maxOverlap = cellParams[3][2]*cellParams[4][2];
    radius = cellParams[4][2]/2.0;
    deathProb = cellParams[5][2]; // fixed death prob. based on CD8+ cells
    divProb = cellParams[11][2]/2.0; // fixed division prob. based on CD8+ cells
    migrationSpeed = 0; // no migration for generalized cell
    influenceRadius = cellParams[9][2];
    rmax = 1.5*radius*2;
    init_time = init_tstamp;
    location_history.push_back(x);
}

void Lymphoid::initialize_cell_from_file(int cell_state, int cell_list_length, double mean_cancer_cell_cycle_length, double std_cancer_cell_cycle_length, RNG& master_rng)
{
    if (cell_state != 12)
    {
        throw std::runtime_error("Lymphoid cell cannot be initialized with cell state = " + std::to_string(cell_state) + ".");
    }
    state = cell_state;
    runtime_index = cell_list_length;
    birthTime = master_rng.uniform(0, 1 / deathProb);
    mother_uniqueID = -1; // A mother_uniqueID of -1 indicates that the cell was part of the initialization of the tumor.
}

void Lymphoid::inherit(std::vector<double> properties)
{
    divProb = properties[0];
    deathProb = properties[1];
}

std::vector<double> Lymphoid::inheritanceProperties()
{
    return {divProb, deathProb};
}

std::array<double, 3> Lymphoid::proliferate(double dt, RNG& master_rng)
{
    if (!canProlif || state == -1)
    {
        return {0,0,0}; // cannot proliferate because suppressed or dead!
    }
    return prob_proliferate(dt, master_rng);
}


