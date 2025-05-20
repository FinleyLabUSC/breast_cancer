#include "../inc/Environment.h"
#include "../inc/ModelUtil.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>

void Environment::initializeCellsFromFile(std::string filePathway) {

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

                // Create and initialize the new cell
                if (cell_state==6) {
                    tempTrajectory = getTcellTrajectory();
                }

                Cell newCell = Cell({x,y}, cellParams,cell_types[cell_state],tempTrajectory);
                newCell.state = cell_state;
                newCell.runtime_index = cell_list.size();

                if (newCell.state == 3) {
                    std::uniform_real_distribution<double> cellCyclePos(0,newCell.cellCycleLength);
                    newCell.cellCyclePos = cellCyclePos(mt);
                    std::uniform_real_distribution<double> pdl1_distribution(0,1);
                    newCell.pdl1 = pdl1_distribution(mt);
                }
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

        tumorSize();
        save(0, 0);
        countPops_updateTimeSeries();
        recordPopulation(0.0);

        std::cout<<"Model initialized. Populations recorded. "<<std::endl;
    } else {
        std::cout<<"\033[31mCouldn't open mIHC file: "<< filePathway << ". Check file pathway or file name!\033[0m"<<std::endl;
    }
}




void Environment::initializeInVitro() {
    std::uniform_real_distribution<double> cell_location(-1500, 1500);

    int idx = 0;
    for (int i = 0; i < 1000; i++) {
        double x = cell_location(mt);
        double y = cell_location(mt);
        Cell newCell = Cell({x,y}, cellParams,0,tCellPhenotypeTrajectory_1);
        std::uniform_real_distribution<double> cellCyclePos(0,newCell.cellCycleLength);
        newCell.cellCyclePos = cellCyclePos(mt);
        newCell.runtime_index = cell_list.size();
        cell_list.push_back(newCell);
        ++idx;
    }

}


void Environment::initializeTesting() {
    Cell newCell = Cell({-100,-100}, cellParams,0,tCellPhenotypeTrajectory_1);
    newCell.runtime_index = cell_list.size();
    cell_list.push_back(newCell);

    newCell =Cell({-100,0}, cellParams,0,tCellPhenotypeTrajectory_1);
    newCell.runtime_index = cell_list.size();
    cell_list.push_back(newCell);

    newCell =Cell({-100,100}, cellParams,0,tCellPhenotypeTrajectory_1);
    newCell.runtime_index = cell_list.size();
    cell_list.push_back(newCell);

}

void Environment::initializeCells() {
    /*
     * places the initial tumor, which is a cluster of pure cancer cells
     */

    double radiiCells = envParams[0];
    int q = 1;
    cell_list.push_back(Cell({0.0,0.0}, cellParams, 0, tCellPhenotypeTrajectory_1));
    cell_list.back().runtime_index = q;

    for(int i=1; i<radiiCells; ++i){
        double circumfrence = 2*i*cellParams[4][0]*3.1415;
        double nCells = circumfrence/cellParams[4][0];
        for(int j=0; j<nCells; ++j){
            double x = i * cellParams[4][0] * cos(2 * 3.1415 * j / nCells);
            double y = i * cellParams[4][0] * sin(2 * 3.1415 * j / nCells);
            Cell newCell = Cell({x, y},  cellParams, 0, tCellPhenotypeTrajectory_1);
            std::uniform_real_distribution<double> cellCyclePos(0,newCell.cellCycleLength);
            newCell.cellCyclePos = cellCyclePos(mt);
            newCell.runtime_index = cell_list.size();

            cell_list.push_back(newCell);
            q++;
        }
    }
}


std::vector<std::string> Environment::getTcellTrajectory() {
    if(tCellPhenotypeTrajectory.empty() || tCellPhenotypeTrajectory.size()==0 ) {
        populateTrajectories(tCellTrajectoryPathway);
    }

    size_t phenotypeIdx = getRandomNumber(tCellPhenotypeTrajectory.size());
    return get2dvecrow(tCellPhenotypeTrajectory, phenotypeIdx);
}


void Environment::recruitImmuneCells(double tstep,  size_t step_count) {

    // recruitment is scaled by number of cancer cells
    auto day = static_cast<double>(tstep*steps/24.0);
    if(day < recruitmentDelay){return;}

    int numC = cancerTS[steps - recruitmentDelay*24/tstep];

    for(int i=0; i<immuneCellRecRates.size(); ++i){ // From genParams.py: [0] cd8, [1] macrophages, [2] cd4, [3] nk, [4] mdsc

        double recRate = immuneCellRecRates[i]*static_cast<double>(numC);
        immuneCells2rec[i] += tstep*recRate; //QUESTION: why are we multiplying with the timestep size? Why are we INCREMENTING?
        while(immuneCells2rec[i] >= 1){
            std::array<double, 2> recLoc =recruitmentLocation_random();

            //i==0 represents the idx in immuneCellRecRates vector that is associated with initializing a cd8 t cell
            if(i == 0){ //
                std::vector<std::string> trajec_phenotype = getTcellTrajectory();

                cell_list.push_back(Cell(recLoc, cellParams, immuneCellRecTypes[i],
                                    trajec_phenotype, step_count));
                cell_list.back().runtime_index = cell_list.size()-1;
            }
            else{ // Recruit macrophage (i=1) or CD4 (i=2) or NK (i=3) or MDSC (i=4)
                cell_list.push_back(Cell(recLoc,cellParams,immuneCellRecTypes[i], tCellPhenotypeTrajectory_1));
                cell_list.back().runtime_index = cell_list.size()-1;
                // This tCellPhenotypeTrajectory_1 is not used for adding CD4 or macrophage
            }
            
            immuneCells2rec[i] -= 1;
        }
    }

}

/** This function generates and returns a random location within the square defined by the tumorRadius size.
 * The location is used for placing immune cells randomly within the domain.
 */
std::array<double, 2> Environment::recruitmentLocation_random() {
    double recDist = 200;
    std::uniform_real_distribution<double> cell_location(-1*tumorRadius-recDist, tumorRadius+recDist);
    double x = tumorCenter[0] + cell_location(mt);
    double y = tumorCenter[1] + cell_location(mt);
    return {x,y};
}

/** This function returns a location within the domain, that is a specific distance away from the tumor center, in a circular region.
 *
 */
std::array<double, 2> Environment::recruitmentLocation() {
    /*
     * cells enter a random distance away from the tumor radius
     * cells enter at a random angle from the tumor center
     */
    //std::uniform_real_distribution<double> angle(-1.0, 1.0);
    std::normal_distribution<double> angle(0.0,1.0);
    std::array<double, 2> dx = {angle(mt),
                                angle(mt)};
    double norm = sqrt(dx[0]*dx[0] + dx[1]*dx[1]);

    std::uniform_real_distribution<double> loc(0.0, recDist);
    double distance = loc(mt) + tumorRadius;

    return {distance*(dx[0]/norm), distance*(dx[1]/norm)};
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









// genearte a truncated normal
void Environment::truncatedGaussian(double mean, double stddev, double min, double max) {
    // Create the list of samples from the truncated distribution, round the values
    int freq;

    std::vector<double> dist_samples;// = {10,20,30,30,20,10,10,20};
    int numSamples = 10000;//dist_samples.size();
    std::vector<double> unique_data;
    std::vector<int> frequency_count;
    std::vector<double> empirical_cdf;
    std::normal_distribution<double> gauss_dist(mean,stddev);


    while (dist_samples.size() < numSamples) {
        double rnd = gauss_dist(mt);
       if (rnd >= min && rnd <= max) {
            //dist_samples.push_back(std::round(rnd)); // sample the normal distribution and round the value to the nearest hour.
            dist_samples.push_back(rnd); // sample the normal distribution and round the value to the nearest hour.
        }
    }

    std::ofstream myfile;
    myfile.open ("generatedDistribution.csv");

    for (auto val:dist_samples) {
        myfile<< val << "\n";
    }
    myfile.close();

    // Sort the list, extract disinct values into the unique_data vector
    std::vector<double> temp_distSamples = dist_samples;
    std::sort(temp_distSamples.begin(), temp_distSamples.end()); // sorts the vector
    unique_data = temp_distSamples;

    auto lastEl = std::unique(unique_data.begin(),unique_data.end()); // finds the unique elements.
    unique_data.resize(std::distance(unique_data.begin(),lastEl)); // resize the unique data vector according to how many elements are unique

    // determine how many of each value occurs in the array.
    for(double i : unique_data) {
        auto bounds = std::equal_range(temp_distSamples.begin(),temp_distSamples.end(),i);
        freq = bounds.second - temp_distSamples.begin() - (bounds.first - temp_distSamples.begin());
        frequency_count.push_back(freq);
    }

    unique_data.insert(unique_data.begin(),0.0); // add zero
    empirical_cdf.push_back(frequency_count[0]);

    double tempSum;

    // Save the CDF / return the CDF
    for(int i = 1; i<frequency_count.size();++i) {
        tempSum = (empirical_cdf[i-1] + frequency_count[i]);
        empirical_cdf.push_back(tempSum);
    }

    for(int i = 0; i < empirical_cdf.size(); i++) {
        empirical_cdf[i] = empirical_cdf[i]/numSamples;
    }

    empirical_cdf.insert(empirical_cdf.begin(),0.0);
}