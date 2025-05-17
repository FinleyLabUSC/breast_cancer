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
    pdl1WhenExpressed = cellParams[8][0];
    pdl1Shift = cellParams[9][0];
    pdl1_increment = cellParams[10][0];

    migrationSpeed = 200; // TODO: select more appropriate migration speed for cancer cells. Current: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7547847/ Fig 1C, Parental median

    rmax = 1.5*radius*2;

    cellCycleLength = setCellCycleLength(state,type);
    init_time = init_tstamp;

    cellCyclePos = 0;
    currDivTime = -1;
    prevDivTime = -1;
    cellAge = 0;
    birthTime = 0;

    deathTime = -1;

    location_history.push_back(x);
}



void Cell::cancer_dieFromCD8(std::array<double, 2> otherX, double otherRadius, double kp, double dt) {
    /*
     * die from CTL based on a probability
     * contact required
     */
    if(type != 0){return;}

    if(calcDistance(otherX) <= radius+otherRadius){
        std::uniform_real_distribution<double> dis(0.0,1.0); 
        if(dis(mt) < kp){
            state = -1;
        }
    }
}

void Cell::cancer_dieFromNK(std::array<double,2> otherX, double otherRadius, double kp, double dt) {
    if (type !=0) {return;}

    if(calcDistance(otherX) <= radius+otherRadius){
        std::uniform_real_distribution<double> dis(0.0,1.0);
        if(dis(mt) < kp){
            state = -1;
        }
    }
}

void Cell::cancer_gainPDL1(double dt) {
    /*
     * shift pdl1 value based on influence from CTL and Th
     *  ifn-g is shown to increase PD-L1 expression
     */
    if(type != 0 || pdl1 > 0.0){return;}

    // induced by ifn-g secreting cells
    // posInfluence is Th + CD8 + NK
    double posInfluence = 1 - (1 - influences[4])*(1 - influences[6])*(1 - influences[8]);
    std::uniform_real_distribution<double> dis(0.0,1.0);
    if(dis(mt) < posInfluence*pdl1Shift){
        pdl1 += pdl1_increment;
    }
}
