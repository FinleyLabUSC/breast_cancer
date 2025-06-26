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
        cell_list[i].next_death_prob = cell_list[i].deathProb;

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
        cell_list[i].indirectInteractions(tstep, step_count,rng,temporary_rng1,anti_pd1_TS.back(),binding_rate_pd1_drug);
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
        if (cell_list[i].deathProb != cell_list[i].next_death_prob) {
            cell_list[i].deathProb = cell_list[i].next_death_prob;
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

    int count_num_cd8_proliferation = 0;
    for(int i=0; i<numCells; ++i){
        if (cell_list[i].state != -1) {
            cell_list[i].set_cellAge(step_count); // this function figures out the age of the cell.
            cell_list[i].age(tstep, step_count,rng); // this function figures out if the cell is dying because it's reached it's lifespan
            cell_list[i].proliferationState(anti_ctla4_TS.back());
            std::array<double, 3> newLoc = cell_list[i].proliferate(tstep, rng);

            if(newLoc[2] == 1){

                if(cell_list[i].type == 0) {
                    // if cancer cell has divided, reset the mother cell's cellCyclePos
                    cell_list[i].prevDivTime = cell_list[i].currDivTime;
                    cell_list[i].currDivTime = step_count;
                    cell_list[i].cellCyclePos =0;
                    cell_list[i].canProlif = false;
                }

                if (cell_list[i].type ==3) {
                    // increase count if CD8+ cell proliferated. Testing purposes.
                    count_num_cd8_proliferation++;
                }

                cell_list.push_back(Cell({newLoc[0], newLoc[1]},cellParams,cell_list[i].type));
                cell_list.back().birthTime = model_time;
                cell_list.back().runtime_index = cell_list.size()-1;

                cell_list[cell_list.size() - 1].inherit(cell_list[i].inheritanceProperties());
            }
        }
    }
    record_proliferation(step_count,count_num_cd8_proliferation);
}

void Environment::runCells(double tstep, size_t step_count) {
    neighborInfluenceInteractions(tstep, step_count);
    calculateForces(tstep,step_count);
    internalCellFunctions(tstep, step_count);
}



void Environment::anti_pd1_drug(double tstep, double new_dose) {
    double current_Value = anti_pd1_TS.back() * std::exp(-1*tstep * anti_pd1_decay_rate) + new_dose;
    double temp = (current_Value>0) ? current_Value : 0;
    anti_pd1_TS.push_back(temp);
}

void Environment::anti_ctla4_drug(double tstep, double new_dose) {
    double current_drug_Value = anti_ctla4_TS.back() * std::exp(-1*tstep * anti_ctla4_decay_rate) + new_dose;
    double temp = (current_drug_Value>0) ? current_drug_Value : 0;
    anti_ctla4_TS.push_back(temp);
}

void Environment::anti_pd1(double tstep) {
    int count_cells_inhibited=0;
    for (auto& cell : cell_list) {
        if (cell.type == 3 || cell.type == 4) { // CD8 cells and NK cells
            double dt = 1; // this can be used to fine grain the pdl1 calculations
            double steps = 1; // tstep / dt; // Use a multiplier to convert from hours (current) to whatever units needed.
            for (int i = 0; i < steps; ++i) {
                cell.pd1_drug_bound = anti_pd1_TS.back()/(anti_pd1_TS.back() + binding_rate_pd1_drug);
                cell.pd1_available = 1- cell.pd1_drug_bound;
                count_cells_inhibited++;
            }
        }
    }
}


void  Environment::treatment(int tx_flag) {}


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

