#include "../inc/Environment.h"

void Environment::printStep(double time) {
    int numM = 0;
    int numT8 = 0;
    int numT4_th = 0;
    int numT4_treg = 0;
    int numC = 0;
    int numNK = 0;
    int numMDSC = 0;

    for(auto &cell : cell_list){

        if(cell.type == 0){
            numC++;
        } else if(cell.type == 1){
            numM++;
        } else if(cell.type == 2){
            if (cell.state == 4) {
                numT4_th++;
            } else {
                numT4_treg++;
            }
        } else if(cell.type == 3){
            numT8++;
        } else if (cell.type == 4) {
            numNK++;
        } else if (cell.type == 5) {
            numMDSC++;
        }
    }


    std::cout << std::fixed << std::setprecision(5);
    std::cout << "Time: " << std::setw(10) << (time / 24) << " | cancer: " << std::setw(10) << numC
          << " | cd8: " << std::setw(10) << numT8 << " | cd4 Th: " << std::setw(10) << numT4_th  << " | cd4 Treg: " << std::setw(10) << numT4_treg
          << " | macrophage: " << std::setw(10) << numM <<  " | nk: " << std::setw(10) << numNK << " | mdsc: " << std::setw(10) << numMDSC << std::endl;

}

void Environment::updateTimeSeries() {
    int numT8 = 0;
    int numT4_th = 0;
    int numT4_treg = 0;
    int numC = 0;
    int numNK = 0;
    int numMDSC = 0;

    int m0 = 0;
    int m1 = 0;
    int m2 = 0;

    // TODO update once I've decided on how to handle the phenotypes of the CD8's.  N -> Pro-M -> exh
    for(auto &cell : cell_list){
        if(cell.type == 3){
            numT8++;
        } else if(cell.type == 2) {
            if (cell.state == 4) {
                numT4_th++;
            } else {
                numT4_treg++;
            }
        } else if(cell.type == 0){
            numC++;
        } else if (cell.type == 4) {
            numNK++;
        } else if (cell.type == 5) {
            numMDSC++;
        } else if (cell.type==1) {
            if (cell.state == 0) {
                if (cell.state == 0) { m0++; }
                if (cell.state == 1) { m1++; }
                if (cell.state == 2) { m2++; }
            }
        }
    }

    cancerTS.push_back(numC);
    cd8TS.push_back(numT8);
    cd4_th_TS.push_back(numT4_th);
    cd4_treg_TS.push_back(numT4_treg);
    nkTS.push_back(numNK);
    mdscTS.push_back(numMDSC);
    m0TS.push_back(m0);
    m1TS.push_back(m1);
    m2TS.push_back(m2);

    radiusTS.push_back(tumorRadius);

}
