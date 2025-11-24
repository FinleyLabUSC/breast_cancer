#include "../inc/RS_Cell.h"

long unsigned RS_Cell::cell_counter = 0; // Defines the cell counter obj

// Initialization
RS_Cell::RS_Cell(std::array<double, 2> loc, std::vector<std::vector<double>>& cellParams, int cellType, size_t init_tstamp): unique_cell_ID(cell_counter++)
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
void RS_Cell::initialize_cell_from_file(int state, int cell_list_length, double mean_cancer_cell_cycle_length, double std_cancer_cell_cycle_length, RNG& master_rng)
{
    throw std::runtime_error("You cannot initialize an RS_Cell object from file.");
}

void RS_Cell::migrate_NN(double dt, RNG& master_rng, std::mt19937& temporary_rng)
{
    // The default behavior here is to perform a biased random walk toward the nearest cancer cell
    // Suppressed cells cannot move as well
    if (state == -1 || immuneSynapseFormed || state == 7 || state == 9)
    {
        location_history.push_back(x); // still save the location history for plotting purposes
        return;
    }

    // Define the random unit vector
    double temp_x = master_rng.normal(0,1,temporary_rng);
    double temp_y = master_rng.normal(0,1,temporary_rng);
    std::array<double, 2> dx_random = {temp_x, temp_y};
    dx_random = unitVector(dx_random);
    std::array<double, 2> dx_movement = {0,0};

    // Find the nearest cancer cell
    std::array<double, 2> nearestCancer = {0.0, 0.0};
    bool nearestCancerFound = false;
    double min_distance = 100 * rmax;
    for (auto& otherCell : cancer_neighbors) {
        double dist = calcDistance(otherCell);
        if (dist < min_distance) {
            min_distance = dist;
            nearestCancer = otherCell;
            nearestCancerFound = true;
        }
    }

    if (nearestCancerFound) {
        std::array<double,2> targetCellDirection =  {nearestCancer[0] - x[0], nearestCancer[1] - x[1]};
        if (std::isnan(targetCellDirection[0]) || std::isnan(targetCellDirection[1])) {
            throw std::runtime_error("NaN encountered in targetCellDirection");
        }
        targetCellDirection = unitVector(targetCellDirection);
        if (std::isnan(targetCellDirection[0]) || std::isnan(targetCellDirection[1])) {
            throw std::runtime_error("NaN encountered after normalizing targetCellDirection");
        }

        for(int i = 0; i < 2; ++i){
            dx_movement[i] = migrationBias * targetCellDirection[i] + (1- migrationBias) * dx_random[i];
        }
        dx_movement = unitVector(dx_movement);
    } else {
        // only take random into account
        dx_movement = dx_random;
    }

    // Migrate!
    for(int i=0; i<2; ++i){
        x[i] += dt*migrationSpeed*dx_movement[i];
        if(std::isnan(x[i])){
            throw std::runtime_error("migration NaN");
        }
    }
    location_history.push_back(x);
}

void RS_Cell::indirectInteractions(double tstep, size_t step_count, RNG& master_rng, std::mt19937& temporary_rng, double anti_pd1_concentration, double binding_rate_pd1_drug)
{
    // indirectInteractions must be specified by cell type
}

void RS_Cell::directInteractions(int interactingState, std::array<double, 2> interactingX, std::vector<double> interactionProperties, double tstep, RNG& master_gen, std::mt19937& temporary_rng)
{
    // directInteractions must be specified by cell type
}

std::vector<double> RS_Cell::directInteractionProperties(int interactingState, size_t step_count)
{
    std::cout << "A cell has tried to call directInteractionProperties from RS_Cell!" << std::endl;
    return {}; // default return nothing for direct interactions as only some cells can do this
}

void RS_Cell::differentiate(double dt, RNG& master_rng, std::mt19937& temporary_rng)
{
    // Differentiation mechanics must be specified per cell type
}

void RS_Cell::express_PDL1(double dt)
{

    if (state == 2 || state == 3 || state == 5 || state == 10) // M2, cancer, Treg, MDSC
    {
        double posInfluence = 1 - (1-influences[4])*(1 - influences[6])*(1 - influences[8]);

        if (posInfluence >= threshold_for_pdl1_induction) {
            double pdl1_increase_amount = (posInfluence - threshold_for_pdl1_induction) * pdl1_induction_rate * dt;
            pdl1_expression_level += pdl1_increase_amount;
            pdl1_expression_level = (pdl1_expression_level < max_pdl1_level) ? pdl1_expression_level : max_pdl1_level;
        } else {
            double pdl1_decrease_amount = pdl1_decay * dt;
            pdl1_expression_level -= pdl1_decrease_amount;
            pdl1_expression_level = (pdl1_expression_level > 0) ? pdl1_expression_level : 0;
        }
    }
}

void RS_Cell::express_PD1(double dt, double anti_pd1_concentration, double binding_rate_pd1_drug)
{
    if (state == 6 || state == 8) // CD8s and NK cells
    {
        // interaction with pro-tumor and tumor cells induce PD1. These also happen to be cells that can express PDL1, but that isn't considered here.
        double influence = 1 - (1-influences[2]) * (1-influences[3])*(1-influences[5])*(1-influences[10]); // interactions with m2, cancer, treg, mdsc
        if (influence >= threshold_for_pd1_induction ) {
            double pd1_increase_amount = (influence - threshold_for_pd1_induction) * pd1_induction_rate * dt;
            pd1_expression_level += pd1_increase_amount;
            pd1_expression_level = (pd1_expression_level < max_pd1_level)?pd1_expression_level : max_pd1_level;
        } else {
            double pd1_decrease_amount = pd1_decay_rate * dt;
            pd1_expression_level -= pd1_decrease_amount;
            pd1_expression_level = (pd1_expression_level > 0) ? pd1_expression_level : 0;
        }

        // determines if the cell is "sensitive" to the drug based on how exhausted it is.
        double fraction_pd1_bound_by_drug = sensitivity_to_antiPD1(anti_pd1_concentration,binding_rate_pd1_drug);
        pd1_drug_bound = pd1_expression_level * fraction_pd1_bound_by_drug;
        pd1_available = pd1_expression_level * (1-fraction_pd1_bound_by_drug);
    }
}

void RS_Cell::pdl1_inhibition(std::array<double, 2> otherX, double otherRadius, double otherPDL1, double dt, RNG& master_rng, std::mt19937& temporary_rng)
{
    // NK and CD8 inhibition by PD1:PDL1 occurs the same way via direct contact w/ PDL1-bearing cells
    // This is not an override right now because we consider this behavior "standard" for killing cells
    if (state == 6 || state == 8){
        double distance = calcDistance(otherX);
        if (distance <= radius+otherRadius)
        {
            double percent_bound = std::min(otherPDL1, pd1_available) / pd1_expression_level;
            double rnd = master_rng.uniform(0, 1, temporary_rng);
            if (rnd < percent_bound)
            {
                next_killProb = next_killProb * percent_bound * inhibitory_effect_of_binding_PD1_PDL1;
                next_migrationSpeed = next_migrationSpeed * percent_bound * inhibitory_effect_of_binding_PD1_PDL1;
                next_death_prob = next_death_prob * (1+percent_bound) * inhibitory_effect_of_binding_PD1_PDL1;
            }
        }
    }
}

void RS_Cell::update_indirectProperties(size_t step_count)
{
    // Which properties to update needs to be specified in the cell type
}

void RS_Cell::contact_die(int killer_state, std::array<double, 2> otherX, double otherRadius, double kill_prob, double dt, RNG& master_rng, std::mt19937& temporary_rng)
{
    // If a cell experiences contact death it must be specified in the cell type
}

std::array<double, 3> RS_Cell::proliferate(double dt, RNG& master_rng)
{
    return {0, 0, 0}; // Default return that the cell does not proliferate
}

void RS_Cell::proliferationState(double anti_ctla4_concentration)
{
    // If the cell does not have its own canProlif logic then it cannot proliferate at all
    canProlif = false;
}

void RS_Cell::mutate(RNG& master_rng) {
    // If a cell can mutate it must be specified in the cell type
}

void RS_Cell::inherit(std::vector<double> properties)
{
    // If a cell can create daughter cells, declare the inherit function in the cell type
}

std::vector<double> RS_Cell::inheritanceProperties()
{
    // Default behavior is to return nothing unless specified in the cell type
    return {};
}


// NEW SHARED PROLIFERATION FUNCTIONS

std::array<double, 3> RS_Cell::cycle_proliferate(double dt, RNG& master_rng)
{
    // These cells set canProlif to false under normal conditions
    // When they reach the end of their clock, they divide
    std::array<double, 2> dx = {master_rng.normal(0, 1), master_rng.normal(0, 1)};
    double norm = calcNorm(dx);
    return{
        radius * (dx[0]/norm) + x[0],
        radius * (dx[1]/norm) + x[1],
        1};
}

std::array<double, 3> RS_Cell::prob_proliferate(double dt, RNG& master_rng)
{
    // These cells set canProlif to false if suppressed
    // They can divide at all other times they can divide
    if (master_rng.uniform(0, 1) < divProb)
    {
        std::array<double, 2> dx = {master_rng.normal(0, 1), master_rng.normal(0, 1)};
        double norm = calcNorm(dx);
        return{
            radius * (dx[0]/norm) + x[0],
            radius * (dx[1]/norm) + x[1],
            1};
    }

    return {0, 0, 0};
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
     * CELLS DO NOT DIE WHEN DETERMINING SPATIAL EQUILIBRIUM
     */
    double rand = master_rng.uniform(0,1);
    if(rand < deathProb){
        if (state == 3) {
            // death_type = 0;
        }

        // state = -1;
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
    if (4 * killProb < kill_prob_base && 4 * migrationSpeed < migration_speed_base  && deathProb > 4 * death_prob_base)
    { // anti-PD1 isn't effective when the CD8 T cell is severely suppressed. "Severe" here is relative killProb < 0.25, relative migration < 0.25 and relative death > 4
        return 0;
    }
    return anti_pd1_concentration / (anti_pd1_concentration + binding_rate_pd1_drug);
}

double RS_Cell::sensitivity_to_antiCTLA4() {
    if (4 * killProb < kill_prob_base && 4 * migrationSpeed < migration_speed_base  && deathProb > 4 * death_prob_base)
    { // anti-PD1 isn't effective when the CD8 T cell is severely suppressed. "Severe" here is relative killProb < 0.25, relative migration < 0.25 and relative death > 4
        return 0;
    }
    return 1;
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
    std::string day_dir = saveDir + "/location_history/";
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
    // calculate influence using an exponential decay based on distance from cell center
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

double RS_Cell::Hill_function(double concentration, double EC50, double n) {
    double effect = 1.0 / (1.0 + pow(EC50 / concentration,n));
    return effect;
}