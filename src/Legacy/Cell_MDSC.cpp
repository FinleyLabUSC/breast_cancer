//
// Created by Rebecca Bekker on 12/30/24.
//

#include "../../inc/Cell.h"

void Cell::initialize_MDSC_Cell(std::vector<std::vector<double> > &cellParams, size_t init_tstamp) {
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

}

void Cell::mdsc_gainPDL1(double dt) {
    if (type != 5){return;}

    // induced by ifn-g secreting cells
    // posInfluence is Th + CD8 + NK
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


