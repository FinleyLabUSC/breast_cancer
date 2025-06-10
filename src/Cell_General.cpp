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
Cell::Cell(std::array<double, 2> loc, std::vector<std::vector<double>> &cellParams, int cellType, std::vector<std::string> tCellPhenotypeTrajectory, size_t init_tstamp): unique_cell_ID(cell_counter++) {
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
    pdl1Shift = 0;
    influenceRadius = 0;
    pdl1 = 0;
    pdl1WhenExpressed = 0;
    pdl1_increment = 0;
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

    // BASELINE CHEMOTHERAPY SENSITIVITY
    chemoDamage = 0; // for all cells
    chemoRepair = 0.0002; // for all cells
    chemoAccumulated = 0; // for all cells
    chemoUptake = 0.085;//for all cells

    // Values for immune cells, cancer cells can acquire tolerance to chemo so they have different values.
    chemoTimeThresh = 0;

  //  std::normal_distribution<double> toleranceDistribution{40,5};
    //chemoTolerance = toleranceDistribution(mt); // this needs to be changed to the correct rng.
    chemo_Accumulated_Threshold = 0;
    chemoTolRate = 0;

    antiPDL1_bindingRate = 1;

    if(cellType == 0){
        initialize_Cancer_Cell(cellParams,0);
        mutationProbability_chemo = 0.0001;
        mutationProbability_inherent = 0;
        chemoTimeThresh = 0;
        chemo_Accumulated_Threshold = 0.0000001; // threshold of "low" dose
        chemoTolRate=0.00; // rate at which cancer cells develop tolerance

    } else if(cellType == 1){
        initialize_Macrophage_Cell(cellParams);
    } else if(cellType == 2){
        initialize_CD4_Cell(cellParams);
    } else if(cellType == 3){
        initialize_CD8_Cell(cellParams, tCellPhenotypeTrajectory, init_tstamp);
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
    if(!canProlif){return {0,0,0};}

    if (type == 0) { // CANCER CELL PROLIFERATION
        // place daughter cell a random angle away from the mother cell at a distance
        std::array<double, 2> dx = {master_rng.normal(0,1),master_rng.normal(0,1)};

        double norm = calcNorm(dx);
        return{radius*(dx[0]/norm)+x[0],
               radius*(dx[1]/norm)+x[1],
               1};

    } else if (type == 3 || type == 8) { // CD8 T CELL PROLIFERATION AND NK PROLIFERATION
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

    //CD8 cells only 
    if(type == 3){
        size_t step_alive = step_count -  init_time; 
                    
        if(step_alive < 0){
            std::cout << "ERROR in OPENMP REGION: NEGATIVE LIFESPAN" << std::endl; 
        }

        //either it is in one of 7 end states or still progressing through state of gene evolution 
        if((step_alive*pTypeStateTransition-1) < t_cell_phenotype_Trajectory.size()){

            std::string phenotype = t_cell_phenotype_Trajectory[step_alive*pTypeStateTransition - 1];
            char phenotype_char = phenotype[0]; 
            
            switch(phenotype_char){
            case 'N': 
                if(rand < deathProb){
                    state = -1;
                }
                break; 
            case 'M': 
                if(rand < (deathProb/2)){
                    state = -1;
                }
                break; 
            default: //case 'E' these are cells that are exhausted but haven't been supressed research showing exhausted t cells kill at lower rate     
                if(rand < deathProb){
                    state = -1;
                }
            }
        }
        else {

            char phenotype_char; 
            if (t_cell_phenotype_Trajectory.empty() || (t_cell_phenotype_Trajectory.size() == 0)){
                std::cerr << "WARNING age: t_cell_phenotype_Trajectory is empty!" << std::endl;
                // suppress bad access with Exhausted phenotype
                phenotype_char = 'E'; 
            }
            else{
                phenotype_char = t_cell_phenotype_Trajectory.back()[0];         
            }
            //char phenotype_char = 'E'; //assume all cells out of their timecourse are exhausted 
            switch(phenotype_char){
            case 'N': 
                if(rand < deathProb){
                    state = -1;
                }
                break; 
            case 'M': 
                if(rand < deathProb/2){
                    state = -1;
                }
                break; 
            default: //case 'E' these are cells that are exhausted but haven't been supressed research showing exhausted t cells kill at lower rate
                if(rand < deathProb){
                    state = -1;
                }
            }
        }
    }
    if(rand < deathProb){
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

void Cell::proliferationState() {
    /*
     * cancer cells and CD8 can proliferate
     */
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
        double posInfluence = 1 - (1 - influences[4])*(1 - influences[8]);
        double negInfluence = 1 - (1 - influences[2])*(1 - influences[5])*(1 - influences[10]);

        double scale = posInfluence - negInfluence;
        divProb = scale*divProb_base;
    } else if (type == 8) {
        canProlif = !compressed; // TODO update depending on whether any other cells affect nk cell proliferation
    }else {
        canProlif = false;
    }
}


void Cell::inherit(std::vector<double> properties) { // QUESTION: cancer cells inherit PDL1 expression, do CD8 T cells inherit phenotype???
    /*
     * daughter cells inherit properties from mother
     * even though this is applicable only for cancer cells, I left it this way for ease of running cell functions
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
        pdl1 = properties[0];
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
        return {pdl1, cellCycleLength};

    } else if (state == 4){
        // CD4 helper
        return {};
    } else if (state == 5){
        // CD4 regulatory
        return {};
    } else if (state == 6){
        // CD8 active
        return {};
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

void Cell::indirectInteractions(double tstep, size_t step_count,RNG& master_rng, std::mt19937& temporary_rng) {

    /*
     * after determining total influences on the cell, run the indirect interaction functions
     */
    if (state == 3){
        // cancer
        cancer_gainPDL1(tstep,master_rng, temporary_rng);
        return;
    }
    if (state == 6){
        // CD8 active
        cd8_setKillProb(step_count);
        return;
    }
    if (state == 8) {
        // NK active
        nk_setKillProb(step_count);
        return;
    }
}

void Cell::directInteractions(int interactingState, std::array<double, 2> interactingX, std::vector<double> interactionProperties, double tstep, RNG& master_rng, std::mt19937& temporary_rng) {

    /*
     * when the cell touches another cell, run direct interactions
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
        if(interactingState == 6){
            // interactionProperties = {radius, killProb}
            cancer_dieFromCD8(interactingX, interactionProperties[0], interactionProperties[1], tstep, master_rng, temporary_rng);
        }
        if (interactingState == 8){
            cancer_dieFromNK(interactingX, interactionProperties[0], interactionProperties[1], tstep, master_rng, temporary_rng);
        }
        return;
    } else if (state == 4){
        // CD4 helper
        return;
    } else if (state == 5){
        // CD4 regulatory
        return;
    } else if (state == 6){
        // CD8 active
        if(interactingState == 2 || interactingState == 3 || interactingState == 5 ||  interactingState == 10){
            // interactionProperties = {radius, pdl1}
            cd8_pdl1Inhibition(interactingX, interactionProperties[0], interactionProperties[1], tstep, master_rng, temporary_rng);
        }
        return;
    } else if (state == 7){
        // CD8 suppressed
        return;
    } else if (state == 8) {
        // NK active
        // M2 cells, Cancer cells, Tregs and MDSCs can induce PDL1 expression on NK cells.
        if (interactingState == 2 || interactingState == 3 || interactingState == 5 ||  interactingState == 10) {
            nk_pdl1Inhibition(interactingX, interactionProperties[0], interactionProperties[1], tstep, master_rng, temporary_rng);
        }
        return;
    } else if (state == 9) {
        // NK suppressed
        return;
    } else if (state == 10) {
        // MDSC
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
            return {radius, pdl1};
        }
        return {};
    } else if (state == 3){
        // cancer
        if(interactingState == 6 || interactingState ==8){
            return {radius, pdl1};
        }
        return {};
    } else if (state == 4){
        // CD4 helper
        return {};
    } else if (state == 5){
        // CD4 regulatory
        if(interactingState == 6 || interactingState ==8){
            return {radius, pdl1};
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
            return {radius, pdl1};
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
        int selectPropertyToMutate = master_rng.uniform_int(1,1);
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
                // TODO update this to be more granular / fine grained.
                if (master_rng.uniform_int(0,1) == 1) { // PDL1 expression changes
                    if (pdl1 == pdl1WhenExpressed) { // if PDL1 is expressed, down regulate
                        pdl1 = 0;
                    } else { // if PDL1 is not expressed, induce expression.
                        pdl1 = pdl1WhenExpressed;
                    }
                }
                break;
            }

            default: std::cout<<"Mutate Error: trying to mutate property not on the list."<<std::endl;
        }
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

void Cell::printLocations(int cellNum) {
    std::ofstream outFile;
    outFile.open("/Users/rebeccabekker/Documents/Finley_Lab/Project1/Code/Mutation/inSilico/migrationMutation/control/migrationInfo_" + std::to_string(cellNum) + ".csv", std::ios::out);

    for (int i = 0; i < location_history.size();i++) {
        outFile<<location_history[i][0]<<","<<location_history[i][1]<<std::endl;
    }
    outFile.close();
}