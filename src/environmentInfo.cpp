#include "../inc/Environment.h"

void Environment::printStep(double time) {

    std::cout << std::fixed << std::setprecision(5);
    std::cout << "Time: " << std::setw(10) << (time / 24) << " | cancer: " << std::setw(10) << cancerTS.back()
            << " | cd8: " << std::setw(10) << cd8TS.back() << " | cd4 Th: " << std::setw(10) << cd4_th_TS.back()  << " | cd4 Treg: " << std::setw(10) << cd4_treg_TS.back()
            << " | macrophage: " << std::setw(10) << m0TS.back() + m1TS.back() + m2TS.back() <<  " | nk: " << std::setw(10) << nkTS.back() << " | mdsc: " << std::setw(10) << mdscTS.back()
            << " | myeloid: " << std::setw(10) << myeloidTS.back() << " | lymphoid: " << std::setw(10) << lymphoidTS.back() << " | stromal: " << std::setw(10) << stromalTS.back() << std::endl;

}

void Environment::countPops_updateTimeSeries() {
    int numT8 = 0;
    int numT4_th = 0;
    int numT4_treg = 0;
    int numC = 0;
    int numNK = 0;
    int numMDSC = 0;

    int m0 = 0;
    int m1 = 0;
    int m2 = 0;

    int myeloid = 0;
    int lymphoid = 0;
    int stromal = 0;

    for(auto &cell : cell_list){
        if(cell->type == 0){
            numC++;
        } else if (cell->type==1) {
            if (cell->state == 0) { m0++; }
            if (cell->state == 1) { m1++; }
            if (cell->state == 2) { m2++; }
        } else if(cell->type == 2) {
            if (cell->state == 4) {
                numT4_th++;
            } else if (cell->state == 5){
                numT4_treg++;
            }
        } else if (cell->type == 3){
            numT8++;
        } else if (cell->type == 4){
            numNK++;
        } else if (cell->type == 5){
            numMDSC++;
        } else if (cell->type == 6){
            myeloid++;
        } else if (cell->type == 7){
            lymphoid++;
        } else if (cell->type == 8){
            stromal++;
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
    lymphoidTS.push_back(lymphoid);
    myeloidTS.push_back(myeloid);
    stromalTS.push_back(stromal);
    radiusTS.push_back(tumorRadius);
}
