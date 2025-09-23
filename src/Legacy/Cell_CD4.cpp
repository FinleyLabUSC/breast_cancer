#include "../../inc/Cell.h"

void Cell::initialize_CD4_Cell(std::vector<std::vector<double> > &cellParams, size_t init_tstamp) {
    state = 4;

    mu = cellParams[0][1];
    kc = cellParams[1][1];
    damping = cellParams[2][1];
    maxOverlap = cellParams[3][1]*cellParams[4][1];
    radius = cellParams[4][1]/2.0;
    deathProb = cellParams[5][1];
    migrationSpeed = cellParams[6][1];
    kTr = cellParams[7][1];
    influenceRadius = cellParams[8][1];
    migrationBias = cellParams[9][1];

    rmax = 1.5*radius*2;

}

void Cell::cd4_differentiation(double dt, RNG& master_rng, std::mt19937& temporary_gen) {
    if(state != 4){return;}

    double rnd = master_rng.uniform(0,1,temporary_gen);
    // negInfuence is M2 + alive cancer + MDSC
    double negInfluence = 1 - (1 - influences[2])*(1 - influences[3])*(1 - influences[10]);
    if(rnd < kTr*negInfluence){
        state = 5;
        cd4_pdl1_expression_level(dt);
    }
}

void Cell::cd4_pdl1_expression_level(double dt) {
    if (state != 5) {return;}

    // ifn gamma producers (Th, CD8, NK) nare assumed to induce expression of PDL1
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
