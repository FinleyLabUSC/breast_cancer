#include "../inc/RS_Cell.h"

long unsigned RS_Cell::cell_counter = 0; // Defines the cell counter obj

// Initialization
RS_Cell::RS_Cell(std::array<double, 2> loc, std::vector<std::vector<double>>& cellParams, int cellType, size_t init_tstamp)
{
    // This constructor sets all parameters shared between cells
    type = cellType;
    state = 0;
    x = loc;

    radius = 0;
    compressed = false;
    currentOverlap = 0;
    divProb = 0;
    divProb_base = 0;
    deathProb = 0;
    canProlif = false;
    mu = 0;
    kc = 0;
    damping = 0;
    maxOverlap = 0;
    rmax = 0;
    currentForces = {0,0};
    migrationSpeed = 0;
    migrationBias = 0;

    influenceRadius = 0;
    immuneSynapseFormed = false; // this indicates whether a cancer cell and a CD8 or NK are in contact. It affects how the migration and forces are handled.

    for(int i=0; i<influences.size(); ++i){
        influences[i] = 0;
    }
    for(int i=0; i<chemotaxVals.size(); ++i){ // TODO: can get rid of chemotaxis??
        chemotaxVals[i] = 0;
    }

    kTr = 0;
    kM1 = 0;
    kM2 = 0;
    plasticity = 0;
    killProb = 0;
    baseKillProb = 0;
    infScale = 0;

    // for influence distance, assume a soft-cutoff where p(distance) = probTh
    probTh = 0.01;
    pd1_drug_bound = 0;
    pd1_available = 0;
    pd1_expression_level = 0;
    pdl1_expression_level = 0;
    max_pd1_level = 5; // arbitrary
    threshold_for_pd1_induction = 0.08; // arbitrary
    pd1_decay_rate = 0.005; // arbitrary
    pd1_induction_rate = 0.05;
    threshold_for_pdl1_induction = 0.08;
    pdl1_induction_rate = 0.05;
    pdl1_decay = 0.005; // arbitrary
    max_pdl1_level = 5; // arbitrary
    inhibitory_effect_of_binding_PD1_PDL1 = 0.5; // moderate suppression
    mutationProbability_inherent = 0.0;
    mutation_prob_PDL1 = 0.0;
}

// OVERRIDE FUNCTIONS
void RS_Cell::initialize_cell(std::vector<std::vector<double>> &cellParams, size_t init_tstamp)
{
    throw std::runtime_error("You cannot initialize an RS_Cell object.");
}

void RS_Cell::initialize_cell_from_file(int state, int run_time_index, double mean_cancer_cell_cycle_length, double std_cancer_cell_cycle_length, RNG& master_rng)
{

}

void RS_Cell::migrate_NN(double dt, RNG& master_rng, std::mt19937& temporary_rng)
{

}

void RS_Cell::indirectInteractions(double tstep, size_t step_count, RNG& master_rng, std::mt19937& temporary_rng, double anti_pd1_concentration, double binding_rate_pd1_drug)
{

}

void RS_Cell::directInteractions(int interactingState, std::array<double, 2> interactingX, std::vector<double> interactionProperties, double tstep, RNG& master_gen, std::mt19937& temporary_rng)
{

}

std::vector<double> RS_Cell::directInteractionProperties(int interactingState, size_t step_count)
{
    return {0};
}

void RS_Cell::express_PD1()
{

}

void RS_Cell::express_PD1L()
{

}

void RS_Cell::pdl1_inhibition()
{

}

void RS_Cell::update_indirectProperties()
{

}

void RS_Cell::contact_die()
{

}

std::array<double, 3> RS_Cell::proliferate(double dt, RNG& master_rng)
{
    return {0, 0, 0}; // Default return that the cell does not proliferate
}

void RS_Cell::mutate(RNG& master_rng) {

}

// SHARED FUNCTIONS
// FORCE FUNCTIONS
std::array<double, 2> RS_Cell::attractiveForce(std::array<double, 2> dx, double otherRadius) {
    double dxNorm = calcNorm(dx);
    std::array<double, 2> dxUnit = {dx[0]/dxNorm, dx[1]/dxNorm};
    double sij = radius + otherRadius;

    double scaleFactor = mu*(dxNorm - sij)*exp(-kc*(dxNorm - sij)/sij);
    double F0 = dxUnit[0]*scaleFactor;
    double F1 = dxUnit[1]*scaleFactor;

    return {F0, F1};
}

std::array<double, 2> RS_Cell::repulsiveForce(std::array<double, 2> dx, double otherRadius) {
    double dxNorm = calcNorm(dx);
    std::array<double, 2> dxUnit = {dx[0]/dxNorm, dx[1]/dxNorm};
    double sij = radius + otherRadius;

    double scaleFactor = mu*sij*log10(1 + (dxNorm - sij)/sij);

    if (scaleFactor < -maxRepulsiveForce) {
        scaleFactor = -maxRepulsiveForce;
    }

    double F0 = dxUnit[0]*scaleFactor;
    double F1 = dxUnit[1]*scaleFactor;

    return {F0, F1};
}

void RS_Cell::calculateForces(std::array<double, 2> otherX, double otherRadius, int &otherType) {
    /*
     * assume attractive force only between cancer cells
     */
    double distance = calcDistance(otherX);
    if(distance < rmax){
        std::array<double, 2> dx = {(otherX[0]-x[0]),
                                    (otherX[1]-x[1])};
        if(distance < (radius + otherRadius)){
            std::array<double, 2> force = repulsiveForce(dx, otherRadius);
            currentForces[0] += force[0];
            currentForces[1] += force[1];
        } else if(type == 0 && otherType == 0){ // attraction if both are cancer cells
            std::array<double, 2> force = attractiveForce(dx, otherRadius);
            currentForces[0] += force[0];
            currentForces[1] += force[1];
        }
    }
}

void RS_Cell::resolveForces(double dt, RNG& master_rng, std::mt19937& temporary_rng) {

    // Only resolve forces on cells that have not formed immune synapses.
    if (!immuneSynapseFormed) {
        // resolving the forces between cells.
        x[0] += (dt/damping)*currentForces[0];
        x[1] += (dt/damping)*currentForces[1];
    }
    resetForces(master_rng,temporary_rng);
}

void RS_Cell::resetForces(RNG& master_rng, std::mt19937& temporary_rng) {
    /*
     * resets forces with a slight randomizing factor
     */
    double temp_x = master_rng.uniform(-1,1,temporary_rng);
    double temp_y = master_rng.uniform(-1,1,temporary_rng);
    currentForces = {temp_x,temp_y};
}

void RS_Cell::determine_neighboringCells(std::array<double,2> otherX, int otherCell_runtime_index, int otherCell_state) {
    double dis = calcDistance(otherX);

    if (dis <= 10*rmax) { // check distance between the cells. if they're within the required distance,
        neighbors.push_back(otherCell_runtime_index); // add runtime_index to neighbors

        if (type != 0 && otherCell_state==3) { // if the original cell is immune and the other cell is cancer, add to a different list.
             cancer_neighbors.push_back(otherX);
        }
    }
}

// OVERLAP FUNCTIONS
void RS_Cell::calculateOverlap(std::array<double, 2> otherX, double otherRadius) {
    /*
     * sum up the total overlap distance between the cell and all other cells
     */
    double distance = calcDistance(otherX);
    if(distance < radius + otherRadius){
        currentOverlap += radius + otherRadius - distance;
    }
}

void RS_Cell::resetOverlap() {
    currentOverlap = 0;
}

void RS_Cell::isCompressed() { // :
    compressed = currentOverlap > maxOverlap;

    // This forces cells that have formed immune synapses to not proliferate. The assumption is that these cells (because they're in contact with others) would be affected by contact-inhibition and thus not proliferating.
    if (immuneSynapseFormed) {
        compressed = true;
    }
    resetOverlap();
}

// Cell behavior functions

void RS_Cell::set_cellAge(size_t step_count) {
    cellAge = step_count - birthTime;
}

void RS_Cell::age(double dt, size_t step_count,  RNG& master_rng) {
    /*
     * cells die based on a probability equal to 1/lifespan
     */
    double rand = master_rng.uniform(0,1);
    if(rand < deathProb){
        if (state == 3) {
            death_type = 0;
        }

        state = -1;
    }
}

// CELL INFLUENCE
void RS_Cell::addInfluence(std::array<double, 2> otherX, double otherInfluence, int otherState) {

    /*
     * determine influence based on distance for each cell state
     * total influence is 1
     */
    if(otherState == -1){return;}

    influences[otherState] = 1 - (1 - influences[otherState])*(1 - calcInfDistance(calcDistance(otherX), otherInfluence));
}

void RS_Cell::clearInfluence() {

    for(double & influence : influences){
        influence = 0;
    }
}

void RS_Cell::set_cell_cycle_length(double cell_cycle_length) {
    if (state == -1) {
        std::cout<<"Error: assigning cell cycle length to dead cell!"<<std::endl;
    }
    if (type == 0 && state == 3) {
        cellCycleLength = cell_cycle_length;
    }
}

double RS_Cell::sensitivity_to_antiPD1(double anti_pd1_concentration, double binding_rate_pd1_drug) {
    if (4 * killProb < kill_prob_base && 4 * migrationSpeed < migration_speed_base  && deathProb > 4 * death_prob_base) {
        return 0; // anti-PD1 isn't effective when the CD8 T cell is severely suppressed. "Severe" here is relative killProb < 0.25, relative migration < 0.25 and relative death > 4
    }else {
        return anti_pd1_concentration / (anti_pd1_concentration + binding_rate_pd1_drug);
    }
}

double RS_Cell::sensitivity_to_antiCTLA4() {
    if (4 * killProb < kill_prob_base && 4 * migrationSpeed < migration_speed_base  && deathProb > 4 * death_prob_base) {
        return 0; // anti-PD1 isn't effective when the CD8 T cell is severely suppressed. "Severe" here is relative killProb < 0.25, relative migration < 0.25 and relative death > 4
    }else {
        return 1;
    }
}

void RS_Cell::resetImmuneSynapse() {
    // We assume that NK - Cancer and CD8 - Cancer synapses exist for an hour at most. So they're reset every time step.
    if (state ) {
        immuneSynapseFormed = false;
    }

}

// OTHER FUNCTIONS

void RS_Cell::printLocations(std::string saveDir) {

    std::ofstream myFile;
    std::string day_dir = saveDir + "\\location_history\\";
    std::string str = "mkdir -p " + day_dir;
    const char *command = str.c_str();
    std::system(command);

    myFile.open(day_dir+ std::to_string(unique_cell_ID)+".csv");

    for (auto & pos : location_history) {
        myFile << pos[0] << "," << pos[1] << std::endl;
    }


    myFile.close();
}

double RS_Cell::calcDistance(std::array<double, 2> otherX) {
    double d0 = (otherX[0] - x[0]);
    double d1 = (otherX[1] - x[1]);

    return sqrt(d0*d0 + d1*d1);
}

double RS_Cell::calcInfDistance(double dist, double xth) {
    /*
     * calculate influence using an exponential decay based on distance from cell center
     */
    double alpha = -log2(probTh);
    double lambda = alpha*0.693/xth;

    return exp(-lambda*dist);
}

double RS_Cell::calcNorm(std::array<double, 2> dx){
    return sqrt(dx[0]*dx[0] + dx[1]*dx[1]);
}

std::array<double, 2> RS_Cell::unitVector(std::array<double, 2> v) {
    double norm = calcNorm(v);

    return {v[0]/norm, v[1]/norm};
}

void RS_Cell::updateRunTimeIndex(int index) {
    runtime_index = index;
}