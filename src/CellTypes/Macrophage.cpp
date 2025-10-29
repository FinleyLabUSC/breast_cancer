#include "../inc/RS_Cell.h"

// The child cell class is initialized w/ the parent constructor
// This allows access to specific pre-set variables
Macrophage::Macrophage(std::array<double, 2> loc, std::vector<std::vector<double>>& cellParams, int cellType, size_t init_tstamp): RS_Cell(loc, cellParams, cellType, init_tstamp)
{
    state = 0;
    mu = cellParams[0][3];
    kc = cellParams[1][3];
    damping = cellParams[2][3];
    maxOverlap = cellParams[3][3]*cellParams[4][3];
    radius = cellParams[4][3]/2.0;
    deathProb = cellParams[5][3];
    migrationSpeed = cellParams[6][3];
    kM1 = cellParams[7][3];
    kM2 = cellParams[8][3];
    influenceRadius = cellParams[9][3];
    migrationBias = cellParams[10][3];
    rmax = 1.5*radius*2;
    location_history.push_back(x);
}

void Macrophage::initialize_cell_from_file(int cell_state, int cell_list_length, double mean_cancer_cell_cycle_length, double std_cancer_cell_cycle_length, RNG& master_rng)
{
    if (cell_state != 0 && cell_state != 1 && cell_state != 2)
    {
        throw std::runtime_error("Macrophage cannot be initialized with cell_state = " + std::to_string(cell_state) + ".");
    }
    state = cell_state;
    runtime_index = cell_list_length;
    birthTime = master_rng.uniform(0, 1 / deathProb);
    mother_uniqueID = -1; // A mother_uniqueID of -1 indicates that the cell was part of the initialization of the tumor.

    // M0 and M1 macrophages don't express PDL1; M2s do
    if (state == 2)
    {
        pdl1_expression_level = master_rng.uniform(0,max_pdl1_level);
    }
}

void Macrophage::indirectInteractions(double tstep, size_t step_count, RNG& master_rng, std::mt19937& temporary_rng, double anti_pd1_concentration, double binding_rate_pd1_drug)
{
    // Only M2 macrophages express PDL1
    if (state == 2)
    {
        express_PDL1(tstep);
    }
}

std::vector<double> Macrophage::directInteractionProperties(int interactingState, size_t step_count)
{
    // Only M2 macrophages interact w/ CD8s and NKs
    if (state == 2 && (interactingState == 6 || interactingState == 8))
    {
        return {radius, pdl1_expression_level};
    }
    return {};
}

void Macrophage::differentiate(double dt, RNG& master_rng, std::mt19937& temporary_rng)
{
    // All macrophage types can differentiate from one to another
    double posInfluence = 1 - (1 - influences[4])*(1 - influences[6])*(1 - influences[8]);
    double negInfluence = 1 - (1 - influences[2])*(1 - influences[3])*(1 - influences[5])*(1 - influences[10]);
    double p1 = kM1*posInfluence;
    double p2 = kM2*negInfluence;
    auto p0 = static_cast<double>(state == 0); // p0 is only nonzero for M0s
    double sum = p0 + p1 + p2;

    std::array<double, 3> probs = {p0/sum, (p0 + p1)/sum, (p0 + p1 + p2)/sum};
    int choice = 0;
    double rnd = master_rng.uniform(0,1,temporary_rng);
    for(int i=0; i<3; ++i){
        if(rnd > probs[i]){choice++;}
    }
    state = choice;
    if(state == 1){ // M1 macrophages don't express PDL1
        pdl1_expression_level = 0;
    } else if(state == 2){ // M2 macrophages do
        express_PDL1(dt);
    }
}
