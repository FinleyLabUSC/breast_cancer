#include <omp.h>

#include "../inc/Environment.h"
#include "../inc/ModelUtil.h"


namespace fs = std::filesystem; 

Environment::Environment(std::string folder, std::string set, int base_seed): rng(base_seed) {

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
    dose_anti_pd1 = envParams[5]; // anti PD-1 dose
    dose_anti_ctla4 = envParams[6];// Default should be: 0.6667nM according to calculation in the "Documentation" document.

    tumorCenter = {0,0};
    tumorRadius = 0;
    necroticRadius = 0;

    steps = 0;
    day = 0;

    dt = 0.005;
}


void Environment::simulate(double tstep, int tx, int met, int num_anti_pd1_doses, int num_anti_ctla4_doses, double bind_rate_pd1_drug) {

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

    anti_pd1_treatment_schedule.resize(num_anti_pd1_doses);
    for(int i = 0; i < num_anti_pd1_doses; ++i){
        anti_pd1_treatment_schedule[i]=0.125 + 7.0 * i;
    }

    int anti_pd1_count_num_dose = 0;
    anti_pd1_TS.push_back(0.0);
    bool anti_pd1_on = false;
    anti_pd1_decay_rate = 0.00193;
    binding_rate_pd1_drug = bind_rate_pd1_drug;

    anti_ctl4_treatment_schedule.resize(num_anti_ctla4_doses);
    for(int i = 0; i < num_anti_ctla4_doses; ++i){
        anti_ctl4_treatment_schedule[i]=0.125 + 7.0 * i;
    }

    int anti_ctla4_count_num_dose = 0;
    anti_ctla4_TS.push_back(0.0);
    bool anti_ctla4_on = false;
    anti_ctla4_decay_rate = -log(0.5)/(14.7 * 24.0); // Half  life of approx 14.7 days ref: 10.1111/bcp.12323

    if (tx==2){
        // anti pd1 monotherapy
        anti_pd1_on = true;
    } else if (tx==3) {
        // anti ctla4 monotherapy
        anti_ctla4_on = true;
    } else if (tx==4) {
        // combination
        anti_pd1_on = true;
        anti_ctla4_on = true;
    }

    record_drug((steps * tstep)/24, tx); // saves the drug concentration

   //std::vector<std::string> metLabels = {"../mihc/121223 P9msP61-2 #06 NT2.5 LM Lungs 16d-1_met_2_adjusted.csv", "../mihc/121223 P9msP61-2 #25 NeuN NT2.5 LM 5 wks #1_met_5_adjusted.csv"};
   std::vector<std::string> metLabels = {"../mihc/test.csv"};


     initializeCellsFromFile(metLabels[met-1]);

    std::cout << "starting simulations...\n";

     while(tstep*steps/24 <simulationDuration ) { // simulationDuration
         double timePoint = tstep*steps/24;

         if (!anti_pd1_on) {
             anti_pd1_drug(tstep,0);
         } else {
             if (anti_pd1_count_num_dose < anti_pd1_treatment_schedule.size()) {
                 if (timePoint == anti_pd1_treatment_schedule[anti_pd1_count_num_dose]) {
                     std::cout<<"Dose "<< anti_pd1_count_num_dose + 1<<" of anti PD-1 administered."<<std::endl;
                     anti_pd1_drug(tstep,dose_anti_pd1);
                     anti_pd1(tstep);
                     anti_pd1_count_num_dose++;
                 } else {
                     anti_pd1_drug(tstep,0);
                     anti_pd1(tstep);
                 }
             } else {
                 anti_pd1_drug(tstep,0);
                 anti_pd1(tstep);
             }
         }

         if (!anti_ctla4_on) {
             anti_ctla4_drug(tstep,0);
         } else {
             if (anti_ctla4_count_num_dose<anti_ctl4_treatment_schedule.size()) {
                 if (timePoint == anti_ctl4_treatment_schedule[anti_ctla4_count_num_dose]) {
                     std::cout<<"Dose "<< anti_ctla4_count_num_dose  + 1 <<" anti CTLA-4 administered."<<std::endl;
                     anti_ctla4_drug(tstep,dose_anti_ctla4);
                     anti_ctla4_count_num_dose++;
                 } else {
                     anti_ctla4_drug(tstep,0);
                 }
             } else {
                 anti_ctla4_drug(tstep,0);
             }
         }


        recruitImmuneCells(tstep, tstep*steps);
        runCells(tstep, tstep*steps);
        mutateCells();
        removeDeadCells(); // loops through, removes the dead cells

         shuffleCells(); // shuffles the cells in the list
        updateCell_list(); // loops through, updates the runtimeindex

        tumorSize(); // loops through twice, calculates the tumor center, calculates the furthest distance from the center to a cancer cell.
        steps += 1;
        countPops_updateTimeSeries(); // loops through and saves the count of each cell type
        printStep(steps * tstep); // Prints to screen

        model_time = steps;

        recordPopulation(steps); // saves the count of each cell to a file.
         record_drug((steps * tstep)/24, tx); // saves the drug concentration
        if (fmod(steps * tstep, 12) == 0) { // change the y value depending on how frequently you want the stuff to be saved.
            save(tstep, steps*tstep);
        }

        if (cancerTS.back() == 0) {
            save(tstep, steps*tstep);
            break;
        }
     }
}
