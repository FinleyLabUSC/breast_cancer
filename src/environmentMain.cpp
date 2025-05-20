#include "../inc/Environment.h"
#include "../inc/ModelUtil.h"


namespace fs = std::filesystem; 

void read_trajectory_csv(std::vector<std::vector<std::string>>& csv_data, std::string fpath){
    std::ifstream myfile(fpath);

    //try to read in the csv file
    if(!myfile.is_open()){
        std::cerr << "Error: Unable to open file at " << fpath << std::endl; 
        return; 
    }
    char ch;
    std::string cell;
    std::vector<std::string> current_row;
    bool inside_quotes = false;
    while (myfile.get(ch)) {
        if (ch == '"') {
            inside_quotes = !inside_quotes;
        } else if (ch == ',' && !inside_quotes) {
            current_row.push_back(cell);
            cell.clear();
        } else if (ch == '\n' && !inside_quotes) {
            current_row.push_back(cell);
            cell.clear();
            csv_data.push_back(current_row);
            current_row.clear();
        } else {
            cell += ch;
        }
    }
    // Close the file
    myfile.close(); 
}

std::vector<std::string> get_phenotype(std::vector<std::vector<std::string>>& trajectory_2d_vec){ 
    std::vector<std::string> res; 

    size_t max_col = trajectory_2d_vec[0].size(); 

    for(size_t i = 1; i < trajectory_2d_vec.size(); ++i){ 
        //start at i=1 because the first value is the title, change this if need be
        res.push_back(trajectory_2d_vec[i][max_col - 1]); 
    }

    return res;
}


Environment::Environment(std::string folder, std::string set, std::string tCellTrajectoryPath): mt((std::random_device())()) {

    /*
     * initialize a simulation environment
     * -----------------------------------
     * loads the three parameter files
     *  cell parameters
     *  environment parameters
     *  recruitment parameters
     *
     * set environment variables to their respective values
     */

    saveDir = "./"+folder+"/set_"+set;

    loadParams();

    // cd8RecRate = recParams[0];
    // mRecRate = recParams[1];
    // cd4RecRate = recParams[2];

    for(int i=0; i<3; ++i){
        immuneCellRecRates.push_back(recParams[i]);
        immuneCells2rec.push_back(0.0);
    }

    immuneCellRecTypes = {3, 1, 2, 4, 5}; // {CD8, macrophage, CD4, NK, MDSC} -> same order as RecRates, values for CELL_STATE

    //recDist = recParams[3];
    maxRecCytoConc = recParams[3];
    recruitmentDelay = recParams[4];

    simulationDuration = envParams[1];
    necroticGrowth = envParams[2];
    necroticForce = envParams[3];
    necroticLimit = envParams[4];
    dose_chemo = envParams[5]; // Chemotherapy dose
    dose_ICI = envParams[6];

    tumorCenter = {0,0};
    tumorRadius = 0;
    necroticRadius = 0;

    steps = 0;
    day = 0;

    dt = 0.005;

    tCellTrajectoryPathway = tCellTrajectoryPath;
    // QUESTION: what is happening in these lines? What's the difference between the trajectory files? What's the difference between the trajectory variables?
    // phenotype_out contains files for the WT / control / no treatment case.
    // phenotype_pdl1_out contains files for the treatment case / anti-PDL1 case.
    // attempt to load in t cell trajectory files by iterating over files in tCellTrajectoryPath
    populateTrajectories(tCellTrajectoryPath);

    //attempt to load in t cell trajectory file 
    std::string trajecPath =  "../t_cell_trajectory/1Tcell_Sim_ABM.csv";

    std::vector<std::vector<std::string>> trajec_csv_1;

     read_trajectory_csv(trajec_csv_1, trajecPath);
    // // if we only care about the phenotype we can simply create a std::vector<std::string or char> of phenotypes 
    std::vector<std::string> phenotype_trajec_1 = get_phenotype(trajec_csv_1);
    tCellPhenotypeTrajectory_1 = phenotype_trajec_1;
}

void Environment::populateTrajectories(std::string tCellTrajectoryPath) {
    if(tCellPhenotypeTrajectory.empty() || tCellPhenotypeTrajectory.size()==0 ) {
        printf("Attempting to read trajectories located in %s...\n", tCellTrajectoryPath.c_str());

        for (const auto& entry : std::filesystem::directory_iterator(tCellTrajectoryPath)){
            //make sure we are only reading csv files
            if(entry.path().extension().string() == ".csv"){

                std::vector<std::vector<std::string>> trajec_csv;

                read_trajectory_csv(trajec_csv, entry.path().string());
                std::vector<std::string> ptype = get_phenotype(trajec_csv);
                tCellPhenotypeTrajectory.push_back(ptype);
            }
            else{
                std::cout << entry.path().string() << " cannot be read as a trajectory, if this is unexpected please confirm the integrity of the .csv file \n";
            }
        }
    }
    printf("Finished reading trajectories located in %s...\n", tCellTrajectoryPath.c_str());

}

void Environment::simulate(double tstep, int tx, int met,int numDoses) {
    /*
     * initializes and runs a simulation
     * ---------------------------------
     * place initial tumor
     * run simulation loop
     *  recruit immune cells
     *  run cell functions
     * ends once time limit is reached or there are no more cancer cells
     */
    model_time = 0;
    
    chemo_treatment_schedule.resize(numDoses);
    for(int i = 0; i < numDoses; i++){
        chemo_treatment_schedule[i]=0.125 + 7.0 * i;
    }
    
    int chemo_count_num_dose = 0;
    chemoTS.push_back(0.0);
    bool chemoOn = false;
    chemotherapy_decay_rate = 0.005; // assume chemo has a half life of 24 hours. (in line with Cess 10.1002/cso2.1032 and 10.3389/fphys.2020.00319

    ici_treatment_schedule = {10.625};//{1.00};//, 7.125, 14.125};
    int ici_count_num_dose = 0;
    ICI_TS.push_back(0.0);
    bool ICI_on = false;
    ICI_decay_rate = 0.00193;


    if (tx == 1){
        chemoOn = true;
    } else if (tx==2){
        chemoOn = true;
        ICI_on = true;
    }

    std::vector<std::string> metLabels = {"../mihc/121223 P9msP61-2 #06 NT2.5 LM Lungs 16d-1_met_2_adjusted.csv", "../mihc/121223 P9msP61-2 #25 NeuN NT2.5 LM 5 wks #1_met_5_adjusted.csv"};

    initializeCellsFromFile(metLabels[met-1]);
    


    std::cout << "starting simulations...\n";


    while(tstep*steps/24 <simulationDuration ) { // simulationDuration
        double timePoint = tstep*steps/24;
        recruitImmuneCells(tstep, tstep*steps);

        runCells(tstep, tstep*steps);
        // if chemo is off
        if (!chemoOn) {
            chemotherapy_drug(tstep,0); // update the chemo value
        } else { // if chemo is on
            if (timePoint == chemo_treatment_schedule[chemo_count_num_dose]) {
                std::cout<<"Chemo administered"<<std::endl;
                chemotherapy_drug(tstep,dose_chemo);
                chemotherapy(tstep);
                chemo_count_num_dose++;
            } else {
                chemotherapy(tstep);
                chemotherapy_drug(tstep,0);
            }
        }
        if (!ICI_on) {
            immune_checkpoint_inhibitor_drug(tstep,0);
        } else {
            if (timePoint == ici_treatment_schedule[ici_count_num_dose]) {
                std::cout<<"ICI administered"<<std::endl;
                immune_checkpoint_inhibitor_drug(tstep,dose_ICI);
                immune_checkpoint_inhibitor(tstep);
                ici_count_num_dose++;
            } else {
                    immune_checkpoint_inhibitor_drug(tstep,0);
                    immune_checkpoint_inhibitor(tstep);
            }
        }

        // mutateCells(chemoOn);

        removeDeadCells(); // loops through, removes the dead cells
        updateCell_list(); // loops through, updates the runtimeindex

        tumorSize(); // loops through twice, calculates the tumor center, calculates the furtherest distance from the center to a cancer cell.

        //necrosis(tstep); // Note for future use: if necroses is turned back on & kills cells, confirm the order of the removeDeadCells, updateCell_list, tumorSize and necrosis functions.

        steps += 1;
        countPops_updateTimeSeries(); // loops through and saves the count of each cell type
        printStep(steps * tstep); // Prints to screen

        model_time = steps;

        recordPopulation(steps*tstep); // saves the count of each cell to a file.

        save(tstep, steps*tstep);

        // if (fmod(steps * tstep, 24) == 0) {
        //     // save every simulation day
        //     save(tstep, steps*tstep);
        // }


        if (cancerTS.back() == 0) {
            save(tstep, steps*tstep);
            break;
        }

    }
}
