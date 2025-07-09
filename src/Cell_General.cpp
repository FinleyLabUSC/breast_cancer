#include "../inc/Environment.h"

#include "../inc/Cell.h"

/*
 * CELL TYPES
 * 0 - cancer
 * 1 - macrophage
 * 2 - CD4
 * 3 - CD8
 *
 * CELL STATES
 * -1 - dead
 * 0 - M0
 * 1 - M1
 * 2 - M2
 * 3 - alive (cancer)
 * 4 - Th
 * 5 - Treg
 * 6 - active (CD8)
 * 7 - suppressed (CD8)
 *
 */

long unsigned Cell::cell_counter = 0;  // Define the static variable here

// INITIALIZE CELL TYPE
Cell::Cell(std::array<double, 2> loc, std::vector<std::vector<double>> &cellParams, int cellType, size_t init_tstamp): unique_cell_ID(cell_counter++) {
    /*
     * initialize all parameters to 0
     * set parameters based on cellType
     */

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

    for(int i=0; i<influences.size(); ++i){
        influences[i] = 0;
    }
    for(int i=0; i<chemotaxVals.size(); ++i){
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
    threshold_for_pd1_induction = 0.5; // arbitrary
    pd1_decay_rate = 0.025; // arbitrary

    threshold_for_pdl1_induction = 0.5;
    pdl1_induction_rate = 0.05;
    pdl1_decay = 0.025; // arbitrary
    max_pdl1_level = 5; // arbitrary

    inhibitory_effect_of_binding_PD1_PDL1 = 0.5; // moderate suppression

    if(cellType == 0){
        initialize_Cancer_Cell(cellParams,0);
        mutationProbability_inherent = 0.01;
        mutation_prob_PDL1 = 0.01;

    } else if(cellType == 1){
        initialize_Macrophage_Cell(cellParams);
    } else if(cellType == 2){
        initialize_CD4_Cell(cellParams);
    } else if(cellType == 3){
        initialize_CD8_Cell(cellParams, init_tstamp);
    } else if (cellType ==4 ) {
        initialize_NK_Cell(cellParams);
    } else if (cellType == 5 ) {
        initialize_MDSC_Cell(cellParams);
    }
    else{
        throw std::runtime_error("Cell::Cell -> unavailable cell type");
    }
}


// FORCE FUNCTIONS
// from Osborne 2017
std::array<double, 2> Cell::attractiveForce(std::array<double, 2> dx, double otherRadius) {
    double dxNorm = calcNorm(dx);
    std::array<double, 2> dxUnit = {dx[0]/dxNorm, dx[1]/dxNorm};
    double sij = radius + otherRadius;

    double scaleFactor = mu*(dxNorm - sij)*exp(-kc*(dxNorm - sij)/sij);
    double F0 = dxUnit[0]*scaleFactor;
    double F1 = dxUnit[1]*scaleFactor;

    return {F0, F1};
}

std::array<double, 2> Cell::repulsiveForce(std::array<double, 2> dx, double otherRadius) {
    double dxNorm = calcNorm(dx);
    std::array<double, 2> dxUnit = {dx[0]/dxNorm, dx[1]/dxNorm};
    double sij = radius + otherRadius;

    double scaleFactor = mu*sij*log10(1 + (dxNorm - sij)/sij);
    double F0 = dxUnit[0]*scaleFactor;
    double F1 = dxUnit[1]*scaleFactor;

    return {F0, F1};
}

void Cell::calculateForces(std::array<double, 2> otherX, double otherRadius, int &otherType) {
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

void Cell::resolveForces(double dt, std::array<double, 2> &tumorCenter, double &necroticRadius, double &necroticForce, RNG& master_rng, std::mt19937& temporary_rng) {

    // resolving the forces between cells.
    x[0] += (dt/damping)*currentForces[0];
    x[1] += (dt/damping)*currentForces[1];

    resetForces(master_rng,temporary_rng);
}

void Cell::resetForces(RNG& master_rng, std::mt19937& temporary_rng) {
    /*
     * resets forces with a slight randomizing factor
     */
    double temp_x = master_rng.uniform(-1,1,temporary_rng);
    double temp_y = master_rng.uniform(-1,1,temporary_rng);
    currentForces = {temp_x,temp_y};
}

// void Cell::neighboringCells(std::array<double, 2> otherX, int otherID){
//
//     /*
//      * determine which cells are within 10*maximum interaction distance
//      * stores the index in cell_list (runtime_index of the cell in Environment) of the neighboring cells
//      */
//     double dis = calcDistance(otherX);
//     if(dis <= 10*rmax){
//         neighbors.push_back(otherID);
//     }
// }
//
// void Cell::neighboringCancerCells(std::array<double, 2> otherX, int otherState){
//
//     /*
//      * determine which cancer cells are within 10*maximum interaction distance
//      * stores the index in cell_list (in Environment) of the neighboring cells
//      */
//     double dis = calcDistance(otherX);
//     if (type != 0 && otherState == 3) { // If the cell under consideration is an immune cell, and the other is a cancer cell
//         if(dis <= 10*rmax){
//             cancer_neighbors.push_back(otherX);
//         }
//     }
// }

void Cell::determine_neighboringCells(std::array<double,2> otherX, int otherCell_runtime_index, int otherCell_state) {
    double dis = calcDistance(otherX);

    if (dis <= 10*rmax) { // check distance between the cells. if they're within the required distance,
        neighbors.push_back(otherCell_runtime_index); // add runtime_index to neighbors

        if (type != 0 && otherCell_state==3) { // if the original cell is immune and the other cell is cancer, add to a different list.
             cancer_neighbors.push_back(otherX);
        }
    }
}



// OVERLAP FUNCTIONS
void Cell::calculateOverlap(std::array<double, 2> otherX, double otherRadius) {
    /*
     * sum up the total overlap distance between the cell and all other cells
     */
    double distance = calcDistance(otherX);
    if(distance < radius + otherRadius){
        currentOverlap += radius + otherRadius - distance;
    }
}

void Cell::resetOverlap() {
    currentOverlap = 0;
}

void Cell::isCompressed() { // :
    compressed = currentOverlap > maxOverlap;
    resetOverlap();
}

// CELL BEHAVIOR FUNCTIONS
// UPDATED PROLIFERATION CODE -> CANCER CELLS HAVE CELL CYCLE LENGTHS, NOT DIVISION PROBABILITIES
std::array<double, 3> Cell::proliferate(double dt, RNG& master_rng) {
    // positions 0 and 1 are cell location
    // position 2 is boolean didProliferate?
    if(!canProlif || state==-1){return {0,0,0};}

    if (type == 0) { // CANCER CELL PROLIFERATION
        // place daughter cell a random angle away from the mother cell at a distance
        std::array<double, 2> dx = {master_rng.normal(0,1),master_rng.normal(0,1)};

        double norm = calcNorm(dx);
        return{radius*(dx[0]/norm)+x[0],
               radius*(dx[1]/norm)+x[1],
               1};

    } else if (type == 3) { // CD8 T CELL PROLIFERATION
        if(master_rng.uniform(0,1) < divProb){
            // place daughter cell a random angle away from the mother cell at a distance of radius
            std::array<double, 2> dx = {master_rng.normal(0,1),master_rng.normal(0,1)};
            double norm = calcNorm(dx);
            return{radius*(dx[0]/norm)+x[0],
                   radius*(dx[1]/norm)+x[1],
                   1};
        } else{
            return {0,0,0};
        }
    }
    return {0,0,0};
}


void Cell::set_cellAge(size_t step_count) {
     cellAge = step_count - birthTime;
}

void Cell::age(double dt, size_t step_count,  RNG& master_rng) {
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

void Cell::migrate_NN(double dt, RNG& master_rng, std::mt19937& temporary_rng) {
    // Dead cells, suppressed CD8's, suppressed NKs don't move

    if(state == -1 || state == 7 || state == 9) {return;}
     if (type == 0) { // Cancer
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
    } else {

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

        double dist;

        for (auto& otherCell : cancer_neighbors) {
            dist = calcDistance(otherCell);
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

        for(int i=0; i<2; ++i){
            x[i] += dt*migrationSpeed*dx_movement[i];
            if(std::isnan(x[i])){
                throw std::runtime_error("migration NaN");
            }
        }
    }
}

void Cell::proliferationState(double anti_ctla4_concentration) {
    /*
     * Only Cancer cells and CD8 can proliferate
     */
    if (state==-1) {return;} // dead cells can't proliferate

    if(type == 0){ // Cancer cells
        if (!(state == -1 || compressed)) {
            cellCyclePos++;
        }
        if (static_cast<double>(cellCyclePos) > cellCycleLength && state == 3 && !compressed) {
            canProlif = true;
        }

    } else if(type == 3){ // CD8 T cells
        // CTLs -> presence of Th and NK promotes their proliferation, M2 and Treg decrease it
        // assume CTLs need IL-2 from Th to proliferate
        canProlif = !(state == 7 || compressed);

        divProb = cd8_setProliferationScale(anti_ctla4_concentration)*divProb_base;
    } else {
        canProlif = false;
    }

}


void Cell::inherit(std::vector<double> properties) {
    /*
     * daughter cells inherit properties from mother
     */
    if(state == 0){
        // M0 macrophages
        return;
    } else if (state == 1){
        // M1 macrophages
        return;
    } else if (state == 2){
        // M2 macrophages
        return;
    } else if (state == 3){
        // cancer
        pdl1_expression_level = properties[0];
        cellCycleLength = properties[1];
        return;
    } else if (state == 4){
        // CD4 helper
        return;
    } else if (state == 5){
        // CD4 regulatory
        return;
    } else if (state == 6){
        // CD8 active
        pd1_expression_level = properties[0];

        killProb = properties[1];
        migrationSpeed = properties[2];
        divProb = properties[3];
        deathProb = properties[4];
        return;
    } else if (state == 7){
        // CD8 suppressed
        return;
    } else if (state == 8){
        // NK active
        return;
    } else if (state == 9){
        // NK suppressed
        return;
    } else if (state == 10){
        // MDSC
        return;
    } else {
        std::cout<<"Error: unknown cell state: " << state << std::endl;
    }
}

std::vector<double> Cell::inheritanceProperties() {
    /*
     * returns the properties that go into Cell::inherit
     */
    if(state == 0){
        // M0 macrophages
        return {};
    } else if (state == 1){
        // M1 macrophages
        return {};
    } else if (state == 2){
        // M2 macrophages
        return {};
    } else if (state == 3){
        // cancer
        return {pdl1_expression_level, cellCycleLength};

    } else if (state == 4){
        // CD4 helper
        return {};
    } else if (state == 5){
        // CD4 regulatory
        return {};
    } else if (state == 6){
        // CD8 active
        return {pd1_expression_level, killProb, migrationSpeed, divProb,deathProb};
    } else if (state == 7){
        // CD8 suppressed
        return {};
    } else if (state == 8){
        // NK active
        return {};
    } else if (state == 9){
        // NK suppressed
        return {};
    } else if (state == 10){
        // MDSC
        return {};
    } else {
        std::cout<<"Error: unknown cell state: " << state << std::endl;
    }

    return {};
}

void Cell::indirectInteractions(double tstep, size_t step_count,RNG& master_rng, std::mt19937& temporary_rng, double anti_pd1_concentration, double binding_rate_pd1_drug) {

    /*
     * after determining total influences on the cell, run the indirect interaction functions
     */
    if (state == 3){
        // cancer
        cancer_gainPDL1(tstep);
        return;
    }
    if (state == 2) {
        macrophage_pdl1_expression_level(tstep);
    }
    if (state == 5) {
        cd4_pdl1_expression_level(tstep);
        return;
    }
    if (state == 6){
        // CD8 active
        cd8_pd1_expression_level(tstep,anti_pd1_concentration, binding_rate_pd1_drug);
        cd8_update_properties_indirect();
        return;
    }
    if (state == 8) {
        // NK active
        nk_pd1_expression_level(tstep,anti_pd1_concentration, binding_rate_pd1_drug);
        nk_update_properties_indirect(step_count);
        return;
    }
    if (state == 10) {
        // MDSC cells
        mdsc_gainPDL1(tstep);
        return;
    }
}

void Cell::directInteractions(int interactingState, std::array<double, 2> interactingX, std::vector<double> interactionProperties, double tstep, RNG& master_rng, std::mt19937& temporary_rng) {

    /*
     * when the cell touches another cell, run direct interactions
     */

    if (state == 3){
        // cancer
        if(interactingState == 6){
            // interactionProperties = {radius, killProb}
            cancer_dieFromCD8(interactingX, interactionProperties[0], interactionProperties[1], tstep, master_rng, temporary_rng);
        }
        if (interactingState == 8){
            cancer_dieFromNK(interactingX, interactionProperties[0], interactionProperties[1], tstep, master_rng, temporary_rng);
        }
        return;
    }
    if (state == 6){
        // CD8 active
        if(interactingState == 2 || interactingState == 3 || interactingState == 5 ||  interactingState == 10){
            // interactionProperties = {radius, pdl1}
            cd8_pdl1Inhibition(interactingX, interactionProperties[0], interactionProperties[1], tstep, master_rng, temporary_rng);
        }
        return;
    }
    if (state == 8) {
        // NK active
        // M2 cells, Cancer cells, Tregs and MDSCs can induce PDL1 expression on NK cells.
        if (interactingState == 2 || interactingState == 3 || interactingState == 5 ||  interactingState == 10) {
            nk_pdl1Inhibition(interactingX, interactionProperties[0], interactionProperties[1], tstep, master_rng, temporary_rng);
        }
        return;
    }
}


std::vector<double> Cell::directInteractionProperties(int interactingState, size_t step_count) {

    /*
     * returns the properties that go into Cell::directInteractions
     */
    if(state == 0){
        // M0 macrophages
        return {};
    } else if (state == 1){
        // M1 macrophages
        return {};
    } else if (state == 2){
        // M2 macrophages
        if(interactingState == 6 || interactingState ==8){
            return {radius, pdl1_expression_level};
        }
        return {};
    } else if (state == 3){
        // cancer
        if(interactingState == 6 || interactingState ==8){
            return {radius, pdl1_expression_level};
        }
        return {};
    } else if (state == 4){
        // CD4 helper
        return {};
    } else if (state == 5){
        // CD4 regulatory
        if(interactingState == 6 || interactingState ==8){
            return {radius, pdl1_expression_level};
        }
        return {};
    } else if (state == 6){
        // CD8 active
        if(interactingState == 3){
            return {radius, killProb};
        }
        return{};
    } else if (state == 7){
        // CD8 suppressed
        return {};
    } else if (state == 8) {
        // NK active
        if (interactingState==3){
            return {radius, killProb};
        }
        return {};
    } else if (state == 9) {
        // NK suppressed
        return {};
    } else if (state == 10) {
        // MDSC
        if(interactingState == 6 || interactingState ==8){
            return {radius, pdl1_expression_level};
        }
        return {};
    }
    return {};
}

// DIFFERENTIATION
void Cell::differentiate(double dt, RNG& master_rng, std::mt19937& temporary_rng) {
    /*
     * runs either macrophage differentiation or cd4 differentiation
     */
    if(type == 1){
        macrophage_differentiation(dt, master_rng, temporary_rng);
    }
    if (type == 2) {
        cd4_differentiation(dt,master_rng, temporary_rng);
    }
}

// CELL INFLUENCE
void Cell::addInfluence(std::array<double, 2> otherX, double otherInfluence, int otherState) {

    /*
     * determine influence based on distance for each cell state
     * total influence is 1
     */
    if(otherState == -1){return;}

    influences[otherState] = 1 - (1 - influences[otherState])*(1 - calcInfDistance(calcDistance(otherX), otherInfluence));
}

void Cell::clearInfluence() {

    for(double & influence : influences){
        influence = 0;
    }
}

void Cell::set_cell_cycle_length(double cell_cycle_length) {
    if (state == -1) {
        std::cout<<"Error: assigning cell cycle length to dead cell!"<<std::endl;
    }
    if (type == 0 && state == 3) {
        cellCycleLength = cell_cycle_length;
    }
}


void Cell::mutate(RNG& master_rng) {
    // properties that can mutate: proliferation rate, migration rate, or PDL1 expression,

    double sampleMutation = master_rng.uniform(0,1);

    // Only alive cancer cells will mutate
    if(type == 0 && state == 3 && sampleMutation < mutationProbability_inherent) {

        // Select what to mutate
        int selectPropertyToMutate = master_rng.uniform_int(0,2);
        switch(selectPropertyToMutate) {
            case 0: {
                // proliferation
                // Increase or decrease according to kappa sampled from a uniform distribution [0, std dev of cell cycle].
                // The position within the cell cycle is changed as well.

                float deltaCellCycle = master_rng.uniform(0,2);
                cellCycleLength *= deltaCellCycle;
                break;
            }
            case 1: {
                // migration rate, incremented or decremented by a percentage (-100%, 100%). If migChange < 1 then decrease speed, if migChange > 1 increase speed
                double migrationDelta = master_rng.uniform(0,2);
                migrationSpeed *= migrationDelta;
                break;
            }
            case 2: {
                // PDL1 expression
                double scale = master_rng.uniform(0,2);
                pdl1_expression_level *= scale;
                break;
            }

            default: std::cout<<"Mutate Error: trying to mutate property not on the list."<<std::endl;
        }
    }
}

/**
 * This function is used to initialize specific properties of the cells when the model is initialized using an mIHC slide.
 * Essentially the assumption is that in this case, cells aren't starting "from scratch". There have already been
 * interactions that gave rise to the current state of the tumor. Cancer cells have been cycling, and T cells & NK cells
 * have been suppressed by the tumor. This function assigns those properties and others.
 *
 * @param cell_state specifies what type of cell we're working with (type = 0, state = 3 - Cancer; type = 3 state = 6 - CD8; type = 4 state = 8 - NK).
 * @param cell_list_length specifies the runtime index of the cell
 * @param mean_cancer_cell_cycle_length
 * @param std_cancer_cell_cycle_length
 * @param master_rng used to sample the appropriate distribution
 */
void Cell::initialize_cell_from_file(int cell_state, int cell_list_length, double mean_cancer_cell_cycle_length, double std_cancer_cell_cycle_length, RNG& master_rng) {

    state = cell_state;
    runtime_index = cell_list_length;
    birthTime = master_rng.uniform(0, 1 / deathProb);
    mother_uniqueID = -1; // A mother_uniqueID of -1 indicates that the cell was part of the initialization of the tumor.

    // Cancer. Assign the cell cycle length, and position within the cell cycle.
    if (type == 0 && state == 3) {
        cellCycleLength = master_rng.normal(mean_cancer_cell_cycle_length, std_cancer_cell_cycle_length);
        cellCyclePos = master_rng.uniform(0, cellCycleLength);
        pdl1_expression_level = master_rng.uniform(0,max_pdl1_level);
    } else if (type == 3) { // CD8+ T cells
        pd1_expression_level = master_rng.uniform(0,max_pd1_level);

        // Assign properties such that the cell is exhausted.
        deathProb = master_rng.uniform(4 * death_prob_base, 20 * death_prob_base); // High death rate

        migrationSpeed = master_rng.uniform(0, 0.25 * migration_speed_base); // Low migration speed
        divProb = master_rng.uniform(0, 0.25 * divProb_base); // Low division probability
        killProb = master_rng.uniform(0, 0.25 * kill_prob_base); // Low cyctotoxic effect
    } else if (type == 4) { // Natural Killer cells
        pd1_expression_level = master_rng.uniform(0,max_pd1_level);

        deathProb = master_rng.uniform(4 * death_prob_base, 20 * death_prob_base); // High death rate
        migrationSpeed = master_rng.uniform(0, 0.25 * migration_speed_base);
        killProb = master_rng.uniform(0, 0.25 * kill_prob_base); // Low cyctotoxic effect
    } else if (state == 2 || state == 5 || state == 10){
        pdl1_expression_level = master_rng.uniform(0,max_pdl1_level);
    }
}

double Cell::sensitivity_to_antiPD1(double anti_pd1_concentration, double binding_rate_pd1_drug) {
    if (4 * killProb < kill_prob_base && 4 * migrationSpeed < migration_speed_base  && deathProb > 4 * death_prob_base) {
        return 0; // anti-PD1 isn't effective when the CD8 T cell is severely suppressed. "Severe" here is relative killProb < 0.25, relative migration < 0.25 and relative death > 4
    }else {
        return anti_pd1_concentration / (anti_pd1_concentration + binding_rate_pd1_drug);
    }
}

double Cell::sensitivity_to_antiCTLA4() {
    if (4 * killProb < kill_prob_base && 4 * migrationSpeed < migration_speed_base  && deathProb > 4 * death_prob_base) {
        return 0; // anti-PD1 isn't effective when the CD8 T cell is severely suppressed. "Severe" here is relative killProb < 0.25, relative migration < 0.25 and relative death > 4
    }else {
        return 1;
    }
}

// OTHER FUNCTIONS
double Cell::calcDistance(std::array<double, 2> otherX) {
    double d0 = (otherX[0] - x[0]);
    double d1 = (otherX[1] - x[1]);

    return sqrt(d0*d0 + d1*d1);
}

double Cell::calcInfDistance(double dist, double xth) {
    /*
     * calculate influence using an exponential decay based on distance from cell center
     */
    double alpha = -log2(probTh);
    double lambda = alpha*0.693/xth;

    return exp(-lambda*dist);
}

double Cell::calcNorm(std::array<double, 2> dx){
    return sqrt(dx[0]*dx[0] + dx[1]*dx[1]);
}

std::array<double, 2> Cell::unitVector(std::array<double, 2> v) {
    double norm = calcNorm(v);

    return {v[0]/norm, v[1]/norm};
}

void Cell::updateRunTimeIndex(int index) {
    runtime_index = index;
}
