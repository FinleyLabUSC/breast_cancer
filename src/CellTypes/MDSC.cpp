#include "../inc/RS_Cell.h"

// The child cell class is initialized w/ the parent constructor
// This allows access to specific pre-set variables
MDSC::MDSC(std::array<double, 2> loc, std::vector<std::vector<double>>& cellParams, int cellType, size_t init_tstamp): RS_Cell(loc, cellParams, cellType, init_tstamp)
{
    state = 10;
    mu = cellParams[0][5];
    kc = cellParams[1][5];
    damping = cellParams[2][5];
    maxOverlap = cellParams[3][5]*cellParams[4][5];
    radius = cellParams[4][5]/2.0;
    deathProb = cellParams[5][5];
    migrationSpeed = cellParams[6][5];
    influenceRadius = cellParams[7][5];
    migrationBias = cellParams[8][5];
    rmax = 1.5*radius*2;
    location_history.push_back(x);
}

void MDSC::initialize_cell_from_file(int cell_state, int cell_list_length, double mean_cancer_cell_cycle_length, double std_cancer_cell_cycle_length, RNG& master_rng)
{
    if (cell_state != 10)
    {
        throw std::runtime_error("MDSC cannot be initialized with cell_state =" + std::to_string(cell_state) + ".");
    }
    state = cell_state;
    runtime_index = cell_list_length;
    birthTime = master_rng.uniform(0, 1 / deathProb);
    mother_uniqueID = -1; // A mother_uniqueID of -1 indicates that the cell was part of the initialization of the tumor.
    pdl1_expression_level = master_rng.uniform(0,max_pdl1_level);
}

void MDSC::indirectInteractions(double tstep, size_t step_count, RNG& master_rng, std::mt19937& temporary_rng, double anti_pd1_concentration, double binding_rate_pd1_drug)
{
    // TODO: Update PDL1 expression function to match declaration
    express_PDL1(tstep);
}

std::vector<double> MDSC::directInteractionProperties(int interactingState, size_t step_count)
{
    // MDSCs interact with CD8s and NKs
    if (interactingState == 6 || interactingState ==8)
    {
        return {radius, pdl1_expression_level};
    }
    return {};
}