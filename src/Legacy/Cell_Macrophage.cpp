#include "../../inc/Cell.h"

void Cell::initialize_Macrophage_Cell(std::vector<std::vector<double> > &cellParams, size_t init_tstamp) {
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
}

void Cell::macrophage_differentiation(double dt, RNG& master_rng, std::mt19937& temporary_rng) {
    /*
     * differentiate
     * probability of not differentiating is constant (before scaling)
     *  the closer a macrophage is to other cells, the more likely to differentiate
     * inf-g secreting cells (CTL and Th) promote M1
     * M2, cancer, and Treg promote M2
     */
    if(type != 1){return;}

    // posInfluence is active CD8 + Th + NK
    // negInfluence is M2 + alive cancer + Treg + MDSC
    double posInfluence = 1 - (1 - influences[4])*(1 - influences[6])*(1 - influences[8]);
    double negInfluence = 1 - (1 - influences[2])*(1 - influences[3])*(1 - influences[5])*(1 - influences[10]);
    double p1 = kM1*posInfluence;
    double p2 = kM2*negInfluence;
    auto p0 = static_cast<double>(state == 0);
    double sum = p0 + p1 + p2;

    std::array<double, 3> probs = {p0/sum,
                                   (p0 + p1)/sum,
                                   (p0 + p1 + p2)/sum};
    int choice = 0;
    double rnd = master_rng.uniform(0,1,temporary_rng);
    for(int i=0; i<3; ++i){
        if(rnd > probs[i]){choice++;}
    }
    state = choice;
    if(state == 1){
        pdl1_expression_level = 0;
    } else if(state == 2){
        macrophage_pdl1_expression_level(dt);
    }
}

void Cell::macrophage_pdl1_expression_level(double dt) {

    // CD4, CD8, NK
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
