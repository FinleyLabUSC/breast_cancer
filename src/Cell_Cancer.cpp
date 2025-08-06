#include "../inc/Environment.h"

#include "../inc/Cell.h"

void Cell::initialize_Cancer_Cell(std::vector<std::vector<double>> &cellParams, size_t init_tstamp) {
    state = 3;
    canProlif = false;

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



void Cell::cancer_dieFromCD8(std::array<double, 2> otherX, double otherRadius, double immune_cell_kill_prob, double dt, RNG& master_rng, std::mt19937& temporary_rng) {
    /*
     * die from CTL based on a probability
     * contact required
     */
    if(type != 0){return;}

    if(calcDistance(otherX) <= radius+otherRadius){
        double rnd = master_rng.uniform(0,1,temporary_rng);
        if(rnd < immune_cell_kill_prob){
            next_state = -1;
            death_type = 1;
        }
    }
}

void Cell::cancer_dieFromNK(std::array<double,2> otherX, double otherRadius, double immune_cell_kill_prob, double dt, RNG& master_rng,std::mt19937& temporary_rng) {
    if (type !=0) {return;}

    if(calcDistance(otherX) <= radius+otherRadius){ // this checks whether the interacting cells are in direct contact with each other
        double rnd = master_rng.uniform(0,1,temporary_rng);
        if(rnd < immune_cell_kill_prob){
            next_state = -1;
            death_type = 2;
        }
    }
}

void Cell::cancer_gainPDL1(double dt) {
    /*
     * shift pdl1 value based on influence from CTL and Th and NK cells
     *  ifn-g is shown to increase PD-L1 expression
     */
    if(type != 0){return;}
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
