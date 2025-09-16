#include "../inc/RS_Cell.h"

// The child cell class is initialized w/ the parent constructor
// This allows access to specific pre-set variables
Cancer::Cancer(std::array<double, 2> loc, std::vector<std::vector<double>>& cellParams, int cellType, size_t init_tstamp): RS_Cell(loc, cellParams, cellType, init_tstamp)
{
    // state & proliferation flag
    state = 3;
    canProlif = false;

    // Set parameters
    mu = cellParams[0][0];
    kc = cellParams[1][0];
    damping = cellParams[2][0];
    maxOverlap = cellParams[3][0]*cellParams[4][0];
    radius = cellParams[4][0]/2.0;
    divProb = cellParams[5][0];
    deathProb = cellParams[6][0];
    influenceRadius = cellParams[7][0];
    migrationSpeed = 30; // 30 microns per hour ref 10.1529/biophysj.106.088898. HER2+ breast cancer cells.
    rmax = 1.5*radius*2;
    init_time = init_tstamp;
    cellCyclePos = 0;
    currDivTime = -1;
    prevDivTime = -1;
    cellAge = 0;
    birthTime = 0;
    location_history.push_back(x);
}


void Cancer::initialize_cell_from_file(int cell_state, int cell_list_length, double mean_cancer_cell_cycle_length, double std_cancer_cell_cycle_length, RNG& master_rng)
{
    if (cell_state != 3 && cell_state != -1) // Allow dead cancer cells for now in mIHC initialization
    {
        throw std::runtime_error("Cancer cell cannot be initalized with cell_state = " + std::to_string(cell_state) + ".");
    }
    state = cell_state;
    runtime_index = cell_list_length;
    birthTime = master_rng.uniform(0, 1 / deathProb);
    mother_uniqueID = -1; // A mother_uniqueID of -1 indicates that the cell was part of the initialization of the tumor.

    // Check to make sure the cancer cell is alive
    if (state == 3)
    {
        cellCycleLength = master_rng.normal(mean_cancer_cell_cycle_length, std_cancer_cell_cycle_length);
        cellCyclePos = master_rng.uniform(0, cellCycleLength);
        pdl1_expression_level = master_rng.uniform(0,max_pdl1_level);
    }
}
