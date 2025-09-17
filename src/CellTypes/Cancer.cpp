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

std::array<double, 3> Cancer::proliferate(double dt, RNG& master_rng)
{
    if (!canProlif || state == -1)
    {
        return {0,0,0}; // cannot proliferate because suppressed or dead!
    }

    return cycle_proliferate(dt, master_rng);
}

void Cancer::migrate_NN(double dt, RNG& master_rng, std::mt19937& temporary_rng)
{
    // If the cancer cell is dead or immune synapsed it can't move
    if (state == -1 || immuneSynapseFormed) {return; }
    // do random migration
    double temp_x = master_rng.normal(0,1,temporary_rng);
    double temp_y = master_rng.normal(0,1,temporary_rng);
    std::array<double, 2> rand_unit_vec = unitVector({temp_x,temp_y});
    for(int i=0; i<x.size(); ++i){
        x[i] += dt*migrationSpeed*rand_unit_vec[i];
        if(std::isnan(x[i])){
            throw std::runtime_error("Cancer migration NaN");
        }
    }
    location_history.push_back(x);
}

void Cancer::indirectInteractions(double tstep, size_t step_count, RNG& master_rng, std::mt19937& temporary_rng, double anti_pd1_concentration, double binding_rate_pd1_drug)
{
    express_PDL1(tstep);
}

void Cancer::directInteractions(int interactingState, std::array<double, 2> interactingX, std::vector<double> interactionProperties, double tstep, RNG& master_gen, std::mt19937& temporary_rng)
{
    contact_die(interactingState, interactingX, interactionProperties[0], interactionProperties[1], tstep, master_gen, temporary_rng);
}

void Cancer::contact_die(int killer_state, std::array<double, 2> otherX, double otherRadius, double kill_prob, double dt, RNG& master_rng, std::mt19937& temporary_rng)
{
    // Cells must be in contact to perform contact death
    if (calcDistance(otherX) <= radius + otherRadius)
    {
        int kill_type = 0;
        switch(killer_state)
        {
        case 6: // Death by CD8
            kill_type = 1; break;
        case 8: // Death by NK
            kill_type = 2; break;
        default:
            throw std::runtime_error("Contact killing behavior undefined for cell state " + std::to_string(killer_state) + ".");
        }
        double rnd = master_rng.uniform(0, 1, temporary_rng);
        if (rnd < kill_prob)
        {
            next_state = -1;
            death_type = kill_type;
        }
    }
}

std::vector<double> Cancer::directInteractionProperties(int interactingState, size_t step_count)
{
    // cancer cells can interact w/ CD8s and NK cells
    if (interactingState == 6 || interactingState == 8)
    {
        return {radius, pdl1_expression_level};
    }
    return {};
}

void Cancer::mutate(RNG& master_rng)
{
    if (state == -1){return;} // Only alive cells can mutate

    double sampleMutation = master_rng.uniform(0,1);
    if(sampleMutation < mutationProbability_inherent) {
        // Select what to mutate
        int selectPropertyToMutate = master_rng.uniform_int(0,2);
        switch(selectPropertyToMutate) {
        case 0: { // Proliferation
                // Increase or decrease according to kappa sampled from a uniform distribution [0, std dev of cell cycle].
                // The position within the cell cycle is scale appropriately.
                float deltaCellCycle = master_rng.uniform(0,2);
                cellCycleLength *= deltaCellCycle;
                break;
        }
        case 1: { // PDL1 expression
                double scale = master_rng.uniform(0,2);
                pdl1_expression_level *= scale;
                break;
        }
        default: std::cout<<"Mutate Error: trying to mutate property not on the list."<<std::endl;
        }
    }
}
