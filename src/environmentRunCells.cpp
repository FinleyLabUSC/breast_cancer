#include <omp.h>

#include "../inc/Environment.h"
#include "../inc/ModelUtil.h"
#include <unordered_set>


void Environment::neighborInfluenceInteractions(double tstep, size_t step_count) {

    /*
     * FIRST LOOP
     * - determine neighbors
     * - determine influences on a cell
     * - perform indirect interactions
     *
     * SECOND LOOP
     * - direct interactions
     *
     * THIRD LOOP
     * - differentiate
     */

#pragma omp parallel for schedule(dynamic)
    for(int i=0; i<cell_list.size(); ++i){
        // reset the "next_" members, so that they're synchronized to the correct values for the current time step.
        cell_list[i].next_state = cell_list[i].state;
        cell_list[i].next_killProb = cell_list[i].killProb;
        cell_list[i].next_migrationSpeed = cell_list[i].migrationSpeed;

        // reset neighborhood and influence
        cell_list[i].neighbors.clear();
        cell_list[i].clearInfluence();

        for(int j = 0; j < cell_list.size(); ++j){
            // assume that a cell cannot influence itself
            if(cell_list[i].unique_cell_ID != cell_list[j].unique_cell_ID){
                cell_list[i].determine_neighboringCells(cell_list[j].x,cell_list[j].runtime_index,cell_list[j].state);
                cell_list[i].addInfluence(cell_list[j].x, cell_list[j].influenceRadius, cell_list[j].state);
            }
        }

        unsigned int seed_for_temp_rng1 = rng.get_context_seed(step_count,cell_list[i].unique_cell_ID,1);
        std::mt19937 temporary_rng1(seed_for_temp_rng1);
        cell_list[i].indirectInteractions(tstep, step_count,rng,temporary_rng1);
    }

#pragma omp parallel for schedule(dynamic)
    for(int i=0; i<cell_list.size(); ++i){
        unsigned int seed_for_temp_rng2 = rng.get_context_seed(step_count,cell_list[i].unique_cell_ID,2);
        std::mt19937 temporary_rng2(seed_for_temp_rng2);

        for(auto &c : cell_list[i].neighbors){
            cell_list[i].directInteractions(cell_list[c].state,
                                            cell_list[c].x,
                                            cell_list[c].directInteractionProperties(cell_list[i].state, step_count),
                                            tstep, rng, temporary_rng2);
        }
    }

#pragma omp parallel for schedule(dynamic)
    for(int i=0; i<cell_list.size(); ++i){
        unsigned int seed_for_temp_rng3 = rng.get_context_seed(step_count,cell_list[i].unique_cell_ID,3);
        std::mt19937 temporary_rng3(seed_for_temp_rng3);
        cell_list[i].differentiate(tstep, rng,temporary_rng3);
    }

    // Update the properties that have changed to reflect the new values.
    // Note: this should be updated if any of the other properties of cells are changed in directInteractions
    // or in a way that may cause data race conditions.

    for (int i = 0; i <cell_list.size();++i) {
        if (cell_list[i].state != cell_list[i].next_state) {
            cell_list[i].state = cell_list[i].next_state;
        }
        if (cell_list[i].killProb != cell_list[i].next_killProb) {
            cell_list[i].killProb = cell_list[i].next_killProb;
        }
        if (cell_list[i].migrationSpeed != cell_list[i].next_migrationSpeed) {
            cell_list[i].migrationSpeed = cell_list[i].next_migrationSpeed;
        }
    }
}

void Environment::calculateForces(double tstep, size_t step_count) {
    /*
     * 1. Calculate total force vector for each cell
     * 2. Resolve forces on each cell
     * 3. Determine current overlap for each cell
     * 4. Determine if each cell is compressed
     */

    // divide tstep into smaller steps for solving
    // only solve forces between neighboring cells to improve computation time

    int Nsteps = static_cast<int>(tstep/dt);

    // iterate through Nsteps, calculating and resolving forces between neighbors
    // also includes migration
    for(int q=0; q<Nsteps; ++q) {
        // migrate first
        #pragma omp parallel for schedule(dynamic)
            for(int i=0; i<cell_list.size(); ++i){
                unsigned int seed_for_temp_rng = rng.get_context_seed(step_count,cell_list[i].unique_cell_ID,4);
                std::mt19937 temporary_rng(seed_for_temp_rng);

                cell_list[i].migrate_NN(dt, rng, temporary_rng);
            }

        // calc forces
        #pragma omp parallel for
            for(int i=0; i<cell_list.size(); ++i){
                for(auto &c : cell_list[i].neighbors){
                    cell_list[i].calculateForces(cell_list[c].x, cell_list[c].radius, cell_list[c].type);
                }
            }

        // resolve forces
       #pragma omp parallel for schedule(dynamic)
            for(int i=0; i<cell_list.size(); ++i){
                unsigned int seed_for_temp_rng = rng.get_context_seed(step_count,cell_list[i].unique_cell_ID,5);
                std::mt19937 temporary_rng(seed_for_temp_rng);

                cell_list[i].resolveForces(dt, tumorCenter, necroticRadius, necroticForce,rng, temporary_rng);
            }
    }

        // calculate overlaps and proliferation states
         #pragma omp parallel for
        for(int i=0; i<cell_list.size(); ++i){
            for(auto &c : cell_list[i].neighbors){
                cell_list[i].calculateOverlap(cell_list[c].x, cell_list[c].radius);
            }
            cell_list[i].isCompressed();
        }
}

void Environment::internalCellFunctions(double tstep, size_t step_count) {
    /*
     * cell death via aging
     * cell proliferation
     * remove cell if out of bounds
     */

    int numCells = cell_list.size();

    for(int i=0; i<numCells; ++i){
        if (cell_list[i].state != -1) {
            cell_list[i].set_cellAge(step_count); // this function figures out the age of the cell.
            cell_list[i].age(tstep, step_count,rng); // this function figures out if the cell is dying because it's reached it's lifespan
            cell_list[i].proliferationState();
            std::array<double, 3> newLoc = cell_list[i].proliferate(tstep, rng);

            if(newLoc[2] == 1){
                if(cell_list[i].type == 3){ // CD8 T cells
                    int phenotypeIdx = 0;  // TODO going to change once the trajectory stuff is removed.
                    std::vector<std::string> trajec_phenotype = get2dvecrow(tCellPhenotypeTrajectory, phenotypeIdx);
                    if(trajec_phenotype.empty() || trajec_phenotype.size() == 0){
                        std::cerr << "WARNING INTERNAL CELL FUNCTIONS: t_cell_phenotype_Trajectory is empty!" << std::endl;
                    }
                    cell_list.push_back(Cell({newLoc[0], newLoc[1]},cellParams,cell_list[i].type, trajec_phenotype));
                    cell_list.back().runtime_index = cell_list.size()-1;
                }
                else{ // any other proliferating cell
                    if(cell_list[i].type == 0) {
                        // if cancer cell has divided, reset the mother cell's cellCyclePos
                        cell_list[i].prevDivTime = cell_list[i].currDivTime;
                        cell_list[i].currDivTime = step_count;
                        cell_list[i].cellCyclePos =0;
                        cell_list[i].canProlif = false;
                    }
                    // Note: NK cells don't currently have cell specific type proliferation behaviour.

                    cell_list.push_back(Cell({newLoc[0], newLoc[1]},cellParams,cell_list[i].type, tCellPhenotypeTrajectory_1));
                    cell_list.back().birthTime = model_time;
                    cell_list.back().runtime_index = cell_list.size()-1;
                }
                cell_list[cell_list.size() - 1].inherit(cell_list[i].inheritanceProperties());
            }
        }
    }
}

void Environment::runCells(double tstep, size_t step_count) {
    neighborInfluenceInteractions(tstep, step_count);
    calculateForces(tstep,step_count);
    internalCellFunctions(tstep, step_count);
}

void Environment::chemotherapy_drug(double tstep, double new_dose) {
    double calc = std::exp(-1*tstep * chemotherapy_decay_rate);
    double currentChemoValue = chemoTS.back() * calc + new_dose;
    double temp = (currentChemoValue > 0) ? currentChemoValue : 0;
  //  std::cout<<"chemoTS "<<chemoTS.back()<< " decay " << chemotherapy_decay_rate << " calc " << calc << " new dose "<< new_dose<< std::endl;
    chemoTS.push_back(temp);
}

void Environment::chemotherapy(double tstep) {
    int count_dead = 0;
    for (auto & cell :cell_list) {
        // if the cell is already dead, return
        if (cell.state==-1) {
            continue;
        } else {
            double dt = 1; // this needs to be adjust based on how small a time step we want to do "diffusion"
            double steps = 1;//60*60*tstep/dt;
            for(int i=0; i<steps; ++i) {

                double drugUptake = cell.chemoUptake * chemoTS.back(); // drug uptake
                cell.chemoDamage += dt * (drugUptake - cell.chemoRepair * cell.chemoDamage);

                // Assume only cancer cells acquire tolerance to drug.
                if (cell.type == 3 && cell.state == 0) {
                    cell.chemoAccumulated += dt * drugUptake; // accumulated chemo, used to model tolerance
                    // modeling accumulation of damage, for the induction of tolerance
                    if (cell.chemoAccumulated > cell.chemo_Accumulated_Threshold) {
                        cell.chemoTime+=dt;
                    }
                    // Checking whether tolerance has been induced
                    if (cell.chemoTime > cell.chemoTimeThresh) {
                        cell.chemoTolerance += dt * cell.chemoTolRate;
                    }
                }
            }
            // kill the cell.
            if (cell.chemoDamage > cell.chemoTolerance) {
                cell.state = -1;
                count_dead++;
            }
        }
    }
}

void Environment::immune_checkpoint_inhibitor(double tstep) {
    int count_cells_inhibited = 0;
    for (auto & cell : cell_list) {
        if (cell.state == 2 || cell.state == 3 || cell.state == 5 || cell.state == 10) { // M2, living Cancer, Treg, MDSC all express PDL1
            double dt = 1; // this can be used to fine grain the pdl1 calculations
            double steps = 1; // tstep / dt; // Use a multiplier to convert from hours (current) to whatever units needed.
            for (int i = 0; i < steps; ++i) {
                double bound_PDL1 = cell.antiPDL1_bindingRate * ICI_TS.back();
                cell.pdl1 = (cell.pdl1 - bound_PDL1 > 0) ? cell.pdl1 - bound_PDL1 : 0;
                count_cells_inhibited++;
            }
        }
    }
}

void Environment::immune_checkpoint_inhibitor_drug(double tstep, double new_dose) {
    double current_ICI_Value = ICI_TS.back() * std::exp(-1*tstep * ICI_decay_rate) + new_dose;
    double temp = (current_ICI_Value>0) ? current_ICI_Value : 0;
    ICI_TS.push_back(temp);
}

void  Environment::treatment() {
    // check the time step against the schedules for the drugs.
        // apply the appropriate drug
}


void Environment::mutateCells() {
    for(auto & cell : cell_list) {
        cell.mutate(rng);
    }
}

void Environment::removeDeadCells() {
    // remove dead cells
    std::vector<int> dead;
    for(int i=0; i<cell_list.size(); ++i){
        if(cell_list[i].state == -1){
            dead.push_back(i);
        }
    }
    std::reverse(dead.begin(), dead.end());
    for(auto &i : dead){
        cell_list.erase(cell_list.begin()+i);
    }
}

void Environment::updateCell_list() {
    for(int i=0; i<cell_list.size(); ++i){
        cell_list[i].updateRunTimeIndex(i);
        if(cell_list[i].state == -1){
            throw std::runtime_error("Environment::internalCellFunctions -> dead cell not removed");
        }
    }
}

