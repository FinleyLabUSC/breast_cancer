#include "../inc/Environment.h"
#include "../inc/ModelUtil.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>

void Environment::initializeCellsFromFile(std::string filePathway) {
    RNG initialize_from_file_rng(722); // Use this rng for the initialization. Then the starting point will be exactly the same for all replicates of a condition.
    std::fstream file;
    file.open(filePathway);

    // Check whether file exists
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

                // TODO sample ages for the cells.
                Cell newCell = Cell({x,y}, cellParams,cell_types[cell_state]);
                newCell.initialize_cell_from_file(cell_state,cell_list.size(),mean_cancer_cell_cycle_length, std_cancer_cell_cycle_length,initialize_from_file_rng);

                // Update the cell list
                cell_list.push_back(newCell);
                i++;
            }
        }

        std::cout << "Time: 0 "  << " | cancer: " << std::setw(10) << numCells[3]
        << " | cd8: " << std::setw(10) << numCells[6] << " | cd4: " << std::setw(10) << numCells[4] << " | treg: " << std::setw(10) << numCells[5]  << " | m0: " << std::setw(10) << numCells[0]
        << " | m1: " << std::setw(10) << numCells[1] << " | m2: " << std::setw(10) << numCells[2]  <<  " | nk: " << std::setw(10) << numCells[8] << " | mdsc: " << std::setw(10) << numCells[10] << std::endl;


        std::cout<<"\033[32mModel initialized from mIHC file!\033[0m"<<std::endl;
        file.close(); // Close the input file

        tumorSize(); // always has to be called prior to countPops_updateTimeSeries. This calculates tumorRadius, the other fnx saves tumorRadius.
        save(0, 0);
        countPops_updateTimeSeries();
        recordPopulation(0.0);

        std::cout<<"Model initialized. Populations recorded. "<<std::endl;
    } else {
        std::cout<<"\033[31mCouldn't open mIHC file: "<< filePathway << ". Check file pathway or file name!\033[0m"<<std::endl;
    }
}


void Environment::initializeInVitro() {
    int idx = 0;
    for (int i = 0; i < 30; i++) {
        double x = rng.uniform(-500,500);
        double y = rng.uniform(-500,500);
        Cell newCell = Cell({x,y}, cellParams,0); // TODO update age sampling
        newCell.cellCycleLength = rng.normal(mean_cancer_cell_cycle_length,std_cancer_cell_cycle_length);
        newCell.cellCyclePos = rng.uniform(0,newCell.cellCycleLength);

        newCell.runtime_index = cell_list.size();
        cell_list.push_back(newCell);
        ++idx;
    }

    tumorSize(); // always has to be called prior to countPops_updateTimeSeries. This calculates tumorRadius, the other fnx saves tumorRadius.
    save(0, 0);
    countPops_updateTimeSeries();
    recordPopulation(0.0);
    std::cout<<"Model initialized. Populations recorded. "<<std::endl;

}


void Environment::initializeTesting() {
    Cell newCell = Cell({-100,-100}, cellParams,0);
    newCell.runtime_index = cell_list.size();
    newCell.cellCycleLength = rng.normal(mean_cancer_cell_cycle_length,std_cancer_cell_cycle_length);
    newCell.cellCyclePos = rng.uniform(0,newCell.cellCycleLength);

    cell_list.push_back(newCell);

    newCell =Cell({-100,0}, cellParams,0);
    newCell.runtime_index = cell_list.size();
    newCell.cellCycleLength = rng.normal(mean_cancer_cell_cycle_length,std_cancer_cell_cycle_length);
    newCell.cellCyclePos = rng.uniform(0,newCell.cellCycleLength);
    cell_list.push_back(newCell);

    newCell =Cell({-100,100}, cellParams,0);
    newCell.runtime_index = cell_list.size();
    newCell.cellCycleLength = rng.normal(mean_cancer_cell_cycle_length,std_cancer_cell_cycle_length);
    newCell.cellCyclePos = rng.uniform(0,newCell.cellCycleLength);
    cell_list.push_back(newCell);

}

void Environment::initializeCells() {
    /*
     * places the initial tumor, which is a cluster of pure cancer cells
     */

    double radiiCells = envParams[0];
    int q = 1;
    cell_list.push_back(Cell({0.0,0.0}, cellParams, 0));
    cell_list.back().cellCycleLength = rng.normal(mean_cancer_cell_cycle_length,std_cancer_cell_cycle_length);
    cell_list.back().cellCyclePos = rng.uniform(0,cell_list.back().cellCycleLength);
    cell_list.back().runtime_index = q;

    for(int i=1; i<radiiCells; ++i){
        double circumfrence = 2*i*cellParams[4][0]*3.1415;
        double nCells = circumfrence/cellParams[4][0];
        for(int j=0; j<nCells; ++j){
            double x = i * cellParams[4][0] * cos(2 * 3.1415 * j / nCells);
            double y = i * cellParams[4][0] * sin(2 * 3.1415 * j / nCells);
            Cell newCell = Cell({x, y},  cellParams, 0);
            newCell.cellCycleLength = rng.normal(mean_cancer_cell_cycle_length,std_cancer_cell_cycle_length);
            newCell.cellCyclePos = rng.uniform(0,newCell.cellCycleLength);
            newCell.runtime_index = cell_list.size();

            cell_list.push_back(newCell);
            q++;
        }
    }
}






void Environment::recruitImmuneCells_cancerBirthDeath(double tstep) {
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
            cell_list.push_back(Cell(recLoc,cellParams,immuneCellRecTypes[i]));
            cell_list.back().runtime_index = cell_list.size()-1;

            immuneCells2rec[i] -= 1;
        }
    }

    num_cancer_births = 0;
    num_cancer_deaths = 0;
}




void Environment::recruitImmuneCells_proportionalTumorBurden(double tstep,  size_t step_count) {

    // recruitment is scaled by number of cancer cells
    auto day = static_cast<double>(tstep*steps/24.0);
    if(day < recruitmentDelay){return;}

    int numC = cancerTS[steps - recruitmentDelay*24/tstep];

    for(int i=0; i<immuneCellRecRates.size(); ++i){ // From genParams.py: [0] cd8, [1] macrophages, [2] cd4, [3] nk, [4] mdsc

        double recRate = immuneCellRecRates[i]*static_cast<double>(numC);
        immuneCells2rec[i] += tstep*recRate; //QUESTION: why are we multiplying with the timestep size? Why are we INCREMENTING?
        while(immuneCells2rec[i] >= 1){
            std::array<double, 2> recLoc = generate_random_location_for_immune_recruitment();

             // Recruit CD8 (i==0) or macrophage (i=1) or CD4 (i=2) or NK (i=3) or MDSC (i=4)
                cell_list.push_back(Cell(recLoc,cellParams,immuneCellRecTypes[i]));
                cell_list.back().runtime_index = cell_list.size()-1;

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
        if(c.type == 0) {
            avgX += c.x[0];
            avgY += c.x[1];
            numC += 1;
        }
    }
    avgX /= numC;
    avgY /= numC;

    tumorCenter = {avgX, avgY};

    tumorRadius = 0;
    for(auto& c : cell_list){
        if(c.type == 0){
            tumorRadius = std::max(tumorRadius, c.calcDistance(tumorCenter));
        }
    }
}

void Environment::necrosis(double tstep) {
    int nCancer = 0;
    for(auto &c : cell_list){
        if(c.type == 0){
            ++nCancer;
        }
    }

    necroticRadius += tstep*nCancer*necroticGrowth;
    necroticRadius = std::max(std::min(necroticRadius, tumorRadius-necroticLimit), 0.0);
}

double Environment::calculateDiffusibles(std::array<double, 2> x) {
    double d = 0;
    for(auto & c: cell_list){
        if(c.state == 3){
            // is a cancer cell
            double alpha = -log2(c.probTh);
            double lambda = alpha*0.693/c.influenceRadius;
            double dist = c.calcDistance(x);
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