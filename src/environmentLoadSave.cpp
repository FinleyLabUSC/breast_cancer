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

    myfile.open(day_dir+"/influences.csv");
    for (auto &cell : cell_list) {
        myfile<<cell.unique_cell_ID;
        for (int i = 0 ; i<cell.influences.size();i++) {
            myfile<<","<<cell.influences[i];
        }
        myfile<<std::endl;
    }
    myfile.close();


    myfile.open(day_dir+"/cells.csv");
    for(auto &cell : cell_list){
        //logging cell location, type and state 
        if(cell.type == 3 || cell.type == 4){ //cd8 t cell
            myfile << cell.unique_cell_ID << ","
           << cell.type << ","
           << cell.x[0] << ","
           << cell.x[1] << ","
           << cell.radius << ","
            << cell.state << ","
            << cell.mother_uniqueID << ","
            << cell.pd1_expression_level<< ","
            << cell.pd1_available<< ","
           << cell.migrationSpeed/cell.migration_speed_base << ","
            << cell.killProb/cell.kill_prob_base << ","
            << cell.deathProb/cell.death_prob_base << ","
            << cell.divProb/cell.divProb_base << std::endl;
        }
        else{
            myfile << cell.unique_cell_ID << ","
                << cell.type << ","
                << cell.x[0] << ","
                << cell.x[1] << ","
                << cell.radius << ","
                << cell.state << ","
                << cell.mother_uniqueID << ","
                << cell.pdl1_expression_level <<  ","
            <<cell.cellCycleLength<< ","
              <<-1<< ","
              <<-1<<","
              <<-1<<","
              <<-1<< std::endl;
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

    if (timestamp==0.0) {
        std::string str = "mkdir -p " + saveDir;
        const char *command = str.c_str();
        std::system(command);
        myfile.open(saveDir + "/populations_TS.csv", std::ios_base::app);
        myfile << "Hr" << "," << "m0" << "," << "m1" << "," << "m2" << "," << "c" << "," << "cd4_th" << "," << "cd4_treg" << "," << "cd8" << "," << "nk"<< "," << "mdsc" << ","  << "radius" << std::endl;
        myfile << timestamp << "," << m0TS.back() << "," << m1TS.back() << "," << m2TS.back() << "," <<  cancerTS.back() << "," << cd4_th_TS.back() << "," << cd4_treg_TS.back() << "," << cd8TS.back() << "," << nkTS.back() << "," << mdscTS.back() << ","  << radiusTS.back() << std::endl;
    }

    myfile.open(saveDir + "/populations_TS.csv", std::ios_base::app);
    myfile << timestamp << "," << m0TS.back() << "," << m1TS.back() << "," << m2TS.back() << "," <<  cancerTS.back() << "," << cd4_th_TS.back() << "," << cd4_treg_TS.back() << "," << cd8TS.back() << "," << nkTS.back() << "," << mdscTS.back() << ","  << radiusTS.back() << std::endl;

    myfile.close();
}

void Environment::record_drug(double tstep, int tx_type) {
    if (tx_type == 2 || tx_type == 4) {
        bool file_exists = std::filesystem::exists(saveDir + "/anti_pd1.csv");
        std::ofstream myfile;
        if (!file_exists) {
            std::string str = "mkdir -p " + saveDir;
            const char *command = str.c_str();
            std::system(command);
            myfile.open(saveDir + "/anti_pd1.csv", std::ios_base::app);
            myfile << "Hr" << "," << "anti_pd1"<< std::endl;
            myfile << tstep << "," << anti_pd1_TS.back() << std::endl;
        } else {
            myfile.open(saveDir + "/anti_pd1.csv", std::ios_base::app);
            myfile << tstep << "," << anti_pd1_TS.back() << std::endl;
        }

    }
    if (tx_type == 3 || tx_type == 4) {
        // check if file exists, if not, create it.
        bool file_exists = std::filesystem::exists(saveDir + "/anti_ctla4.csv");
        std::ofstream myfile;
        if (!file_exists) {
            std::string str = "mkdir -p " + saveDir;
            const char *command = str.c_str();
            std::system(command);
            myfile.open(saveDir + "/anti_ctla4.csv", std::ios_base::app);
            myfile << "Hr" << "," << "anti_ctla4"<< std::endl;
            myfile << tstep << "," << anti_ctla4_TS.back() << std::endl;
        } else {
            myfile.open(saveDir + "/anti_ctla4.csv", std::ios_base::app);
            myfile << tstep << "," << anti_ctla4_TS.back() << std::endl;
        }
    }
}

void Environment::record_immuneCount(double tstep, int contact_count) {
    std::ofstream myfile;
    if (tstep==0) {
        std::string str = "mkdir -p " + saveDir;
        const char *command = str.c_str();
        std::system(command);
        myfile.open(saveDir + "/immune_contact.csv", std::ios_base::app);
        myfile << "Hr"  << "," <<  "contact_count"<< std::endl;
    }

    myfile.open(saveDir + "/immune_contact.csv", std::ios_base::app);
    myfile << tstep << "," << contact_count << std::endl;

    myfile.close();
}




void Environment::record_proliferation(double tstep, int prolifCount) {
    std::ofstream myfile;
    if (tstep==0) {
        std::string str = "mkdir -p " + saveDir;
        const char *command = str.c_str();
        std::system(command);
        myfile.open(saveDir + "/cd8_proliferation.csv", std::ios_base::app);
        myfile << "Hr"  << "," <<  "prolif_count"<< std::endl;
    }

    myfile.open(saveDir + "/cd8_proliferation.csv", std::ios_base::app);
    myfile << tstep << "," << prolifCount << std::endl;

    myfile.close();
}


void Environment::record_cancerdeath(double tstep, int count_age_deaths, int count_cd8_contact_deaths, int count_nk_contact_deaths) {
    std::ofstream myfile;
    if (tstep==0) {
        std::string str = "mkdir -p " + saveDir;
        const char *command = str.c_str();
        std::system(command);
        myfile.open(saveDir + "/cancer_death.csv", std::ios_base::app);
        myfile << "Hr" << "," << "age" << "," << "cd8" << "," <<"nk"<< std::endl;
    }

    myfile.open(saveDir + "/cancer_death.csv", std::ios_base::app);
    myfile << tstep << "," << count_age_deaths << "," << count_cd8_contact_deaths  << "," << count_nk_contact_deaths << std::endl;

    myfile.close();
}


