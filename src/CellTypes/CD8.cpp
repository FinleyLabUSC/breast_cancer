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

std::array<double, 3> CD8::proliferate(double dt, RNG& master_rng)
{
    if (!canProlif || state == -1)
    {
        return {0,0,0}; // cannot proliferate because suppressed or dead!
    }

    return prob_proliferate(dt, master_rng);
}

void CD8::indirectInteractions(double tstep, size_t step_count, RNG& master_rng, std::mt19937& temporary_rng, double anti_pd1_concentration, double binding_rate_pd1_drug)
{
    express_PD1(tstep, anti_pd1_concentration, binding_rate_pd1_drug);
    update_indirectProperties(step_count);
}

void CD8::directInteractions(int interactingState, std::array<double, 2> interactingX, std::vector<double> interactionProperties, double tstep, RNG& master_gen, std::mt19937& temporary_rng)
{
    if (interactingState == 2 || interactingState == 3 || interactingState == 5 || interactingState == 10)
    {
        pdl1_inhibition(interactingX, interactionProperties[0], interactionProperties[1], tstep, master_gen, temporary_rng);
    }
}

std::vector<double> CD8::directInteractionProperties(int interactingState, size_t step_count)
{
    // CD8s interact w/ alive cancer cells
    if (interactingState == 3)
    {
        return {radius, killProb};
    }
    return {};
}

void CD8::update_indirectProperties(size_t step_count)
{
    double posInfluence = 1 - (1 - influences[1])*(1 - influences[4]);
    double negInfluence = 1 - (1 - influences[2])*(1 - influences[5])*(1 - influences[10]);
    double scale = posInfluence - negInfluence;

    // TODO: reconsider how these properties are updated
    next_killProb = next_killProb*pow(infScale, scale);
    next_migrationSpeed = next_migrationSpeed*pow(migScale, scale);
}

void CD8::inherit(std::vector<double> properties)
{
    pd1_expression_level = properties[0];
    killProb = properties[1];
    migrationSpeed = properties[2];
    divProb = properties[3];
    deathProb = properties[4];
}

std::vector<double> CD8::inheritanceProperties()
{
    return {pd1_expression_level, killProb, migrationSpeed, divProb, deathProb};
}

void CD8::proliferationState(double anti_ctla4_concentration)
{
    canProlif = !compressed; // CD8 proliferation is only stopped if the cell is compressed

    // scale divProb based on influences & CTLA
    double posInfluence = 1 - (1 - influences[4])*(1 - influences[8]);
    double negInfluence = 1 - (1 - influences[2])*(1 - influences[10]);
    double anti_ctla4_effect = Hill_function(anti_ctla4_concentration,anti_CTLA4_IC50,anti_CTLA4_hill_coeff) * sensitivity_to_antiCTLA4();
    double scale_anti_CTLA4_effect = (divProb > divProb_base) ? 1.1 : 1;
    double effective_antiCTLA4_effect = (anti_ctla4_effect * scale_anti_CTLA4_effect < 1) ? anti_ctla4_effect * scale_anti_CTLA4_effect : 1;
    double ctla_scale = posInfluence * (1 - influences[5]* (1-effective_antiCTLA4_effect)) - negInfluence;
    divProb = ctla_scale * divProb;
}
