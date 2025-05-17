//
// Created by Rebecca Bekker on 12/30/24.
//

#include "../inc/Cell.h"

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
    pdl1WhenExpressed = cellParams[9][5];
    pdl1_increment = cellParams[10][5];
    pdl1 = pdl1WhenExpressed;
    rmax = 1.5*radius*2;

}

void Cell::mdsc_gainPDL1(double dt) {}
