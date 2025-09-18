#include "../inc/Environment.h"
#include "../inc/ModelUtil.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>

#include "RS_Cell.h"

void Environment::initializeCellsFromFile(std::string filePathway) {
    RNG initialize_from_file_rng(722); // Use this rng for the initialization. Then the starting point will be exactly the same for all replicates of a condition.
    std::fstream file;
    file.open(filePathway);

    // Check whether the file exists
    std::vector<int> numCells = {0,0,0,0,0,0,0,0,0,0,0};
    std::vector<int> cell_types = {1,1,1,0,2,2,3,3,4,4,5};

    if (file.is_open()) {
        std::cout<<"\033[32mModel input file opened successfully.\033[0m"<<std::endl;
        std::string line, value;

        // Read in each line
        int i = 0;
        while(getline(file,line)) {
            std::vector<std::string> tempTrajectory;
            std::stringstream ss(line);
            std::vector<std::string> row;

            // Read each value separated by a comma
            while (std::getline(ss, value, ',')) {
                row.push_back(value);
            }

            // Extract each value from the row vector
            double x = std::stod(row[0]); // convert the string to a double value
            double y = std::stod(row[1]);
            int cell_state = stoi(row[2]); // convert the string to an integer value

            if (cell_state!=-10){
                numCells[cell_state]++;

                // TODO: sample ages for the cells.
                addCell({x, y}, cellParams, cell_types[cell_state]);
                // We subtract one from cell_list_length because the cell has ALREADY been added
                cell_list.back()->initialize_cell_from_file(cell_state,cell_list.size() - 1,mean_cancer_cell_cycle_length, std_cancer_cell_cycle_length,initialize_from_file_rng);
                i++;
            }
        }

        std::cout << "Time: 0 "  << " | cancer: " << std::setw(10) << numCells[3]
        << " | cd8: " << std::setw(10) << numCells[6] << " | cd4: " << std::setw(10) << numCells[4] << " | treg: " << std::setw(10) << numCells[5]  << " | m0: " << std::setw(10) << numCells[0]
        << " | m1: " << std::setw(10) << numCells[1] << " | m2: " << std::setw(10) << numCells[2]  <<  " | nk: " << std::setw(10) << numCells[8] << " | mdsc: " << std::setw(10) << numCells[10] << std::endl;


        std::cout<<"\033[32mModel initialized from mIHC file!\033[0m"<<std::endl;
        file.close(); // Close the input file

        tumorSize(); // Always has to be called prior to countPops_updateTimeSeries. This calculates tumorRadius, the other fnx saves tumorRadius.
        save(0, 0);
        countPops_updateTimeSeries();
        recordPopulation(0.0);

        count_cancer_immune_contacts(-1.0);
        std::cout<<"Model initialized. Populations recorded. "<<std::endl;
    } else {
        std::cout<<"\033[31mCouldn't open mIHC file: "<< filePathway << ". Check file pathway or file name!\033[0m"<<std::endl;
    }
}


void Environment::initializeInVitro() {
    int idx = 0;
    for (int i = 0; i < 100; i++) {
        double x = rng.uniform(-20,20);
        double y = rng.uniform(-20,20);
        std::array<double, 2> loc = {x, y};

        // Create cell
        std::shared_ptr<Cancer> newCancer = std::make_shared<Cancer>(loc, cellParams, 0);
        newCancer->cellCycleLength = rng.normal(mean_cancer_cell_cycle_length,std_cancer_cell_cycle_length);
        newCancer->cellCyclePos = rng.uniform(0,newCancer->cellCycleLength);
        newCancer->runtime_index = cell_list.size();
        cell_list.push_back(newCancer);
        ++idx;
    }

    // for (int i = 0; i < 100; i++) {
    //     double x = rng.uniform(-10,10);
    //     double y = rng.uniform(-10,10);
    //     Cell newCell = Cell({x,y}, cellParams,3); // TODO update age sampling
    //
    //     newCell.runtime_index = cell_list.size();
    //     cell_list.push_back(newCell);
    //     ++idx;
    // }

    tumorSize(); // Always has to be called prior to countPops_updateTimeSeries. This calculates tumorRadius, the other fnx saves tumorRadius.
    save(0, 0);
    countPops_updateTimeSeries();
    recordPopulation(0.0);
    count_cancer_immune_contacts(-1.0);
    std::cout<<"Model initialized. Populations recorded. "<<std::endl;

    std::cout << "Time: 0 "  << " | cancer: " << std::setw(10) << cancerTS.back()
       << " | cd8: " << std::setw(10) << cd8TS.back() << " | cd4: " << std::setw(10) << cd4_th_TS.back() << " | treg: " << std::setw(10) << cd4_treg_TS.back()  << " | m0: " << std::setw(10) << m0TS.back()
       << " | m1: " << std::setw(10) << m1TS.back() << " | m2: " << std::setw(10) << m2TS.back()  <<  " | nk: " << std::setw(10) << nkTS.back() << " | mdsc: " << std::setw(10) << mdscTS.back() << std::endl;
}


void Environment::initializeTesting() {
    // A line of cells
    std::array<double, 2> loc1 = {-100, -100};
    std::array<double, 2> loc2 = {-100, 0};
    std::array<double, 2> loc3 = {-100, 100};

    std::shared_ptr<Cancer> newCancer1 = std::make_shared<Cancer>(loc1, cellParams,0);
    newCancer1->runtime_index = cell_list.size();
    newCancer1->cellCycleLength = rng.normal(mean_cancer_cell_cycle_length,std_cancer_cell_cycle_length);
    newCancer1->cellCyclePos = rng.uniform(0,newCancer1->cellCycleLength);
    cell_list.push_back(newCancer1);

    std::shared_ptr<Cancer> newCancer2 = std::make_shared<Cancer>(loc2, cellParams,0);
    newCancer2->runtime_index = cell_list.size();
    newCancer2->cellCycleLength = rng.normal(mean_cancer_cell_cycle_length,std_cancer_cell_cycle_length);
    newCancer2->cellCyclePos = rng.uniform(0,newCancer2->cellCycleLength);
    cell_list.push_back(newCancer2);

    std::shared_ptr<CD8> newCD8 = std::make_shared<CD8>(loc3, cellParams, 0);
    newCD8->runtime_index = cell_list.size();
    newCD8->cellCycleLength = rng.normal(mean_cancer_cell_cycle_length,std_cancer_cell_cycle_length);
    newCD8->cellCyclePos = rng.uniform(0,newCD8->cellCycleLength);
    cell_list.push_back(newCD8);

    tumorSize(); // Always has to be called prior to countPops_updateTimeSeries. This calculates tumorRadius, the other fnx saves tumorRadius.
    save(0, 0);
    countPops_updateTimeSeries();
    recordPopulation(0.0);
    std::cout<<"Model initialized. Populations recorded. "<<std::endl;
    count_cancer_immune_contacts(-1.0);

    std::cout << "Time: 0 "  << " | cancer: " << std::setw(10) << cancerTS.back()
       << " | cd8: " << std::setw(10) << cd8TS.back() << " | cd4: " << std::setw(10) << cd4_th_TS.back() << " | treg: " << std::setw(10) << cd4_treg_TS.back()  << " | m0: " << std::setw(10) << m0TS.back()
       << " | m1: " << std::setw(10) << m1TS.back() << " | m2: " << std::setw(10) << m2TS.back()  <<  " | nk: " << std::setw(10) << nkTS.back() << " | mdsc: " << std::setw(10) << mdscTS.back() << std::endl;

}

void Environment::initializeCells() {
    /*
     * places the initial tumor, which is a cluster of pure cancer cells
     * TODO: Specify these as cancer cells & probably rename the function
     */

    double radiiCells = envParams[0];
    int q = 1;
    std::array<double, 2> loc = {0.0, 0.0};
    std::shared_ptr<Cancer> newCancer = std::make_shared<Cancer>(loc, cellParams, 0);
    cell_list.push_back(newCancer);
    cell_list.back()->cellCycleLength = rng.normal(mean_cancer_cell_cycle_length,std_cancer_cell_cycle_length);
    cell_list.back()->cellCyclePos = rng.uniform(0,cell_list.back()->cellCycleLength);
    cell_list.back()->runtime_index = q;

    for(int i=1; i<radiiCells; ++i){
        double circumfrence = 2*i*cellParams[4][0]*3.1415;
        double nCells = circumfrence/cellParams[4][0];
        for(int j=0; j<nCells; ++j){
            q++;
            double x = i * cellParams[4][0] * cos(2 * 3.1415 * j / nCells);
            double y = i * cellParams[4][0] * sin(2 * 3.1415 * j / nCells);
            std::array<double, 2> loc2 = {x, y};
            std::shared_ptr<Cancer> newCancer2 = std::make_shared<Cancer>(loc2, cellParams, 0);
            newCancer2->cellCycleLength = rng.normal(mean_cancer_cell_cycle_length,std_cancer_cell_cycle_length);
            newCancer2->cellCyclePos = rng.uniform(0,newCancer2->cellCycleLength);
            newCancer2->runtime_index = q;
            cell_list.push_back(newCancer2);
        }
    }

    tumorSize(); // Always has to be called prior to countPops_updateTimeSeries. This calculates tumorRadius, the other fnx saves tumorRadius.
    save(0, 0);
    countPops_updateTimeSeries();
    recordPopulation(0.0);
    std::cout<<"Model initialized. Populations recorded. "<<std::endl;
    count_cancer_immune_contacts(-1.0);

    std::cout << "Time: 0 "  << " | cancer: " << std::setw(10) << cancerTS.back()
       << " | cd8: " << std::setw(10) << cd8TS.back() << " | cd4: " << std::setw(10) << cd4_th_TS.back() << " | treg: " << std::setw(10) << cd4_treg_TS.back()  << " | m0: " << std::setw(10) << m0TS.back()
       << " | m1: " << std::setw(10) << m1TS.back() << " | m2: " << std::setw(10) << m2TS.back()  <<  " | nk: " << std::setw(10) << nkTS.back() << " | mdsc: " << std::setw(10) << mdscTS.back() << std::endl;

}

void Environment::recruitImmuneCells_cancerBirthDeath(double tstep) {
    // TODO: Specify which immune cells are getting recruited when looping!
    // Recruitment of different immune cell types are governed by different cancer-related events.
    // M0, CD4 and MDSC's are recruited proportional to the number of cancer cell "births" in the previous time step.
    // CD8's and NK's are recruited proportional to the number of cancer cell deaths in the previous time step.
    auto day = static_cast<double>(tstep*steps/24.0);
    if(day < recruitmentDelay){return;}

    // CD8 (i==0) or macrophage (i=1) or CD4 (i=2) or NK (i=3) or MDSC (i=4)
    std::vector<int> births_deaths = {num_cancer_deaths, num_cancer_births, num_cancer_births, num_cancer_deaths, num_cancer_births};

    for(int i=0; i<immuneCellRecRates.size(); ++i){ // From genParams.py: [0] cd8, [1] macrophages, [2] cd4, [3] nk, [4] mdsc

        double recRate = immuneCellRecRates[i] * static_cast<double>(births_deaths[i]);
        immuneCells2rec[i] += tstep*recRate;
        while(immuneCells2rec[i] >= 1){
            std::array<double, 2> recLoc = generate_random_location_for_immune_recruitment();

            // Recruit CD8 (i==0) or macrophage (i=1) or CD4 (i=2) or NK (i=3) or MDSC (i=4)
            addCell(recLoc, cellParams, immuneCellRecTypes[i]);
            cell_list.back()->runtime_index = cell_list.size()-1;

            immuneCells2rec[i] -= 1;
        }
    }
    num_cancer_births = 0;
    num_cancer_deaths = 0;
}

void Environment::recruitImmuneCells_proportionalTumorBurden(double tstep,  size_t step_count) {
    // TODO: Specify which immune cells are getting recruited when looping!
    // recruitment is scaled by the number of cancer cells
    auto day = static_cast<double>(tstep*steps/24.0);
    if(day < recruitmentDelay){return;}

    int numC = cancerTS[steps - recruitmentDelay*24/tstep];

    for(int i=0; i<immuneCellRecRates.size(); ++i){ // From genParams.py: [0] cd8, [1] macrophages, [2] cd4, [3] nk, [4] mdsc

        double recRate = immuneCellRecRates[i]*static_cast<double>(numC);
        immuneCells2rec[i] += tstep*recRate; //QUESTION: why are we multiplying with the timestep size? Why are we INCREMENTING?
        while(immuneCells2rec[i] >= 1){
            std::array<double, 2> recLoc = generate_random_location_for_immune_recruitment();

            // Recruit CD8 (i==0) or macrophage (i=1) or CD4 (i=2) or NK (i=3) or MDSC (i=4)
            addCell(recLoc, cellParams, immuneCellRecTypes[i]);
            cell_list.back()->runtime_index = cell_list.size()-1;

            immuneCells2rec[i] -= 1;
        }
    }
}

/** This function generates and returns a random location within the square defined by the tumorRadius + recDist size.
 * The location is used for placing immune cells randomly within the domain.
 */
std::array<double, 2> Environment::generate_random_location_for_immune_recruitment() {
    double recDist = 200;
    double x = tumorCenter[0] + rng.uniform(-1*tumorRadius-recDist, tumorRadius+recDist);
    double y = tumorCenter[1] + rng.uniform(-1*tumorRadius-recDist, tumorRadius+recDist);
    return {x,y};
}

void Environment::tumorSize(){
    tumorCenter = {0,0};
    double avgX = 0;
    double avgY = 0;
    double numC = 0;
    for(auto &c : cell_list){
        if(c->type == 0) {
            avgX += c->x[0];
            avgY += c->x[1];
            numC += 1;
        }
    }
    avgX /= numC;
    avgY /= numC;

    tumorCenter = {avgX, avgY};

    tumorRadius = 0;
    for(auto& c : cell_list){
        if(c->type == 0){
            tumorRadius = std::max(tumorRadius, c->calcDistance(tumorCenter));
        }
    }
}

void Environment::necrosis(double tstep) {
    int nCancer = 0;
    for(auto &c : cell_list){
        if(c->type == 0){
            ++nCancer;
        }
    }

    necroticRadius += tstep*nCancer*necroticGrowth;
    necroticRadius = std::max(std::min(necroticRadius, tumorRadius-necroticLimit), 0.0);
}

double Environment::calculateDiffusibles(std::array<double, 2> x) {
    double d = 0;
    for(auto & c: cell_list){
        if(c->state == 3){
            // is a cancer cell
            double alpha = -log2(c->probTh);
            double lambda = alpha*0.693/c->influenceRadius;
            double dist = c->calcDistance(x);
            d = 1 - (1 - d)*(1 - exp(-lambda*dist));
        }
    }

    return d;
}

/**
 * Shuffles the cells in the cell_list so that there are no spatial artifacts.
 */
void Environment::shuffleCells() {
    std::shuffle(cell_list.begin(), cell_list.end(), rng.getGenerator());
}

void Environment::addCell(std::array<double, 2> loc, std::vector<std::vector<double>>& cellParams, int cellType)
{
    switch (cellType)
    {
    case 0: // Cancer cell
    {
        std::shared_ptr<Cancer> newCancer = std::make_shared<Cancer>(loc, cellParams, cellType);
        cell_list.push_back(newCancer);
        break;
    }
    case 1: // Macrophage
    {
        std::shared_ptr<Macrophage> newMP = std::make_shared<Macrophage>(loc, cellParams, cellType);
        cell_list.push_back(newMP);
        break;
    }
    case 2: // CD4
    {
        std::shared_ptr<CD4> newCD4 = std::make_shared<CD4>(loc, cellParams, cellType);
        cell_list.push_back(newCD4);
        break;
    }
    case 3: // CD8
    {
        std::shared_ptr<CD8> newCD8 = std::make_shared<CD8>(loc, cellParams, cellType);
        cell_list.push_back(newCD8);
        break;
    }
    case 4: // NK
    {
        std::shared_ptr<NK> newNK = std::make_shared<NK>(loc, cellParams, cellType);
        cell_list.push_back(newNK);
        break;
    }
    case 5: // MDSC
    {
        std::shared_ptr<MDSC> newMDSC = std::make_shared<MDSC>(loc, cellParams, cellType);
        cell_list.push_back(newMDSC);
        break;
    }
    default:
        throw std::runtime_error("Cannot create an unknown cell type!");
    }
}