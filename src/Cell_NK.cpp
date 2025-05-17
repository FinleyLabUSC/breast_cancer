#include "../inc/Cell.h"

void Cell::initialize_NK_Cell(std::vector<std::vector<double> > &cellParams, size_t init_tstamp) {
    state = 8;

    mu = cellParams[0][4];
    kc = cellParams[1][4];
    damping = cellParams[2][4];
    maxOverlap = cellParams[3][4]*cellParams[4][4];
    radius = cellParams[4][4]/2.0;
    deathProb = cellParams[5][4];
    migrationSpeed = cellParams[6][4];
    baseKillProb = cellParams[7][4];
    infScale = cellParams[8][4];
    influenceRadius = cellParams[9][4];
    migrationBias = cellParams[10][4];
    divProb_base = cellParams[11][4];

    rmax = 1.5*radius*2;

}

void Cell::nk_pdl1Inhibition(std::array<double, 2> otherX, double otherRadius, double otherpdl1, double dt) {
    // inhibition via direct contact

    if(state != 8){return;}

    double distance = calcDistance(otherX);
    if(distance <= radius+otherRadius){
        std::uniform_real_distribution<double> dis(0.0,1.0);
        if(dis(mt) < otherpdl1){
            state = 9;
            killProb = 0;
            migrationSpeed = 0.0;
        }
    }
}

void Cell::nk_setKillProb(size_t step_count){
    // set kill prob based on influence
    // essentially effects of cytokines

    if(state != 8){return;}

    // posInfluence is M1 + Th
    // negInfluence is M2 + Treg + MDSC
    double posInfluence = 1 - (1 - influences[1])*(1 - influences[4]);
    double negInfluence = 1 - (1 - influences[2])*(1 - influences[5])*(1 - influences[10]);
    double scale = posInfluence - negInfluence;

    killProb = killProb*pow(infScale, scale); //realize the impact of pos and negative influence in coordination with T cell boolean network
}