#include "../inc/Environment.h"

void Environment::loadParams() {
    std::ifstream dataCP(saveDir+"/params/cellParams.csv");
    std::string line;
    while(std::getline(dataCP, line)){
        std::stringstream lineStream(line);
        std::string cell;
        std::vector<double> parsedRow;
        while(std::getline(lineStream, cell, ',')){
            parsedRow.push_back(std::stod(cell));
        }
        cellParams.push_back(parsedRow);
    }
    dataCP.close();

    std::ifstream dataRP(saveDir+"/params/recParams.csv");
    while(std::getline(dataRP, line)){
        std::stringstream lineStream(line);
        std::string cell;
        std::vector<double> parsedRow;
        while(std::getline(lineStream, cell, ',')){
            parsedRow.push_back(std::stod(cell));
        }
        recParams.push_back(parsedRow[0]);
    }
    dataRP.close();

    std::ifstream dataEP(saveDir+"/params/envParams.csv");
    while(std::getline(dataEP, line)){
        std::stringstream lineStream(line);
        std::string cell;
        std::vector<double> parsedRow;
        while(std::getline(lineStream, cell, ',')){
            parsedRow.push_back(std::stod(cell));
        }
        envParams.push_back(parsedRow[0]);
    }
    dataEP.close();
}

void Environment::save(double tstep, double tstamp) {

    std::ofstream myfile;
    std::string day_dir = saveDir + "/cellLists/day_" + std::to_string(day);
    std::string str = "mkdir -p " + day_dir;
    const char *command = str.c_str();
    std::system(command);

    myfile.open(saveDir+"/necroticRadius.csv");
    myfile << necroticRadius << std::endl;
    myfile.close();

    myfile.open(day_dir+"/cells.csv");
    for(auto &cell : cell_list){
        //logging cell location, type and state 
        if(cell.type == 3){ //cd8 t cell
            size_t idx = (tstamp - cell.init_time)*cell.pTypeStateTransition; 
            if(idx > cell.t_cell_phenotype_Trajectory.size() - 1){
                //we are outside of the array and want to get the last 
                std::string phenotype_char; 
                if (cell.t_cell_phenotype_Trajectory.empty() || (cell.t_cell_phenotype_Trajectory.size() == 0)){
                    std::cerr << "WARNING LOGGING: t_cell_phenotype_Trajectory is empty!" << std::endl;
                    phenotype_char = 'E'; 
                }
                else{
                    phenotype_char = cell.t_cell_phenotype_Trajectory.back();         
                }

                myfile << cell.unique_cell_ID << ","
                    << cell.type << ","
                    << cell.x[0] << ","
                    << cell.x[1] << ","
                    << cell.radius << ","
                    << phenotype_char << ","
                    << cell.pdl1 << std::endl;
            }
            else{
                //we are in the array and can index
                std::string pType = cell.t_cell_phenotype_Trajectory[idx]; 
                myfile << cell.unique_cell_ID << ","
                    << cell.type << ","
                    << cell.x[0] << ","
                    << cell.x[1] << ","
                    << cell.radius << ","
                    << pType << ","
                    << cell.pdl1 << std::endl;
            }

        }
        else{
            myfile << cell.unique_cell_ID << ","
                << cell.type << ","
                << cell.x[0] << ","
                << cell.x[1] << ","
                << cell.radius << ","
                << cell.state << ","
                << cell.pdl1 << std::endl;
        }
    }
    myfile.close();

    if(day == 0){
        day++;
        return;}
    ++day;
}

void Environment::recordPopulation(double timestamp) {

    // is there a way that i can check if it's already been created?
    std::ofstream myfile;

    if (timestamp==0) {
        std::string str = "mkdir -p " + saveDir;
        const char *command = str.c_str();
        std::system(command);
        myfile.open(saveDir + "/populations_TS.csv", std::ios_base::app);
        myfile << "Hr" << "," << "m0" << "," << "m1" << "," << "m2" << "," << "c" << "," << "cd4_th" << "," << "cd4_treg" << "," << "cd8" << "," << "nk"<< "," << "mdsc" << ","  << "radius" << std::endl;
    }

    myfile.open(saveDir + "/populations_TS.csv", std::ios_base::app);
    myfile << timestamp << "," << m0TS.back() << "," << m1TS.back() << "," << m2TS.back() << "," <<  cancerTS.back() << "," << cd4_th_TS.back() << "," << cd4_treg_TS.back() << "," << cd8TS.back() << "," << nkTS.back() << "," << mdscTS.back() << ","  << radiusTS.back() << std::endl;

    myfile.close();
}