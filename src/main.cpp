#include <iostream>
#include "../inc/Environment.h"
#include <chrono>

/*
 * CHANGE summing influences to multiplying (1 - inf)
 */

/*
 * MAIN THINGS TO LOOK INTO RIGHT NOW
 * ----------------------------------
 * - impact of macrophages on CD8 infiltration: https://www.pnas.org/doi/10.1073/pnas.1720948115
 * - how M1 and M2 should impact CD8: M2 is easy to find (Petty and Yang, Tumor-associated macrophages: implications in cancer immunotherapy, 2017)
 *      - M2
 *          - reduced CD8 proliferation via cytokine secretion - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6715643/ (media from m2 macrophages)
 *
 *      - M1
 *          - increased T cell infiltration and CD8 activation - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6715643/ (but is this from M1, or simply the removal of M2?)
 *          - conversion to M1 "regulated T-cell response by relieving immunosuppression" - https://pubs.acs.org/doi/full/10.1021/acsami.1c07626
 *          - direct tumor killing - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7751482/
 *   but M1 is a bit more difficult
 *   - changes in CD8 killing effect
 *   - changes in CD8 cytokine secretion
 *   - changes in CD8 proliferation
 *   - changes in CD8 recruitment
 *   - changes in tumor proliferation
 * - M1 cytotoxic effect
 * - CD8 proliferation
 * - cell recruitment
 * - inclusion of Th2
 * - functions of CD4
 *      - Th1
 *          - IFN-y secretion (https://www.nature.com/articles/s41417-020-0183-x)
 *          - IL-2 secretion: drive CD8 effector function, differentiation, and proliferation (https://www.nature.com/articles/s41417-020-0183-x)
 *          - increase number of CD8 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1403291/)
 *          - can replace exogenous IL-2, due to IL-2 production (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1403291/)
 *          - removing Treg is not enough to remove tumor, T help must be provided (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1403291/)
 *          - specific impacts on CD8 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2970736/)
 *              - increased CD8 recruitment to tumor site (IFN-y)
 *              - increased CD8 survival (IL-2)
 *              - increased CD8 proliferation at tumor site (IL-2)
 *              - increased CD8 granzyme B (IL-2)
 *      - Treg
 *          - suppress effects of exogenous IL-2 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1403291/)
 *
 * - suppressed CD8 IFN-y
 *
 * MODEL FRAMEWORK THINGS
 * ----------------------
 * - what is going on with migration bias?
 */

/*
 * MODEL OF TUMOR-IMMUNE INTERACTIONS
 * ----------------------------------
 * - PUT MODEL DESCRIPTION HERE
 */


int main(int argc, char **argv) {
        std::string run_location = argv[1];
        std::string folder = argv[2];
        std::string replicate_number = argv[3];
        std::string pST = argv[4];
        std::string dp_fac = argv[5];
        std::string kp_fac = argv[6];
        int tx = std::stoi(argv[7]);
        int met = std::stoi(argv[8]);
        double binding_rate_pd1_drug = 0.1;
        int repNum = std::stoi(argv[3]);
        double cd8_prolif = std::stod(argv[9]);
        double cd8_death = std::stod(argv[10]);
        double cd8_rec = std::stod(argv[11]);

        // These labels correspond to the type of treatment administered according to the value of tx. If tx = 1, then the simulation uses anti-PD1.
        std::vector<std::string> txLabels = {"control","pd1","ctla4","ici_combo"};
        std::vector<std::string> parameter_levels = {"low","high"};

        int prolif_id = (cd8_prolif==0.08) ? 0 : 1;
        int death_id = (cd8_death==0.005) ? 0 : 1;
        int rec_id = (cd8_rec==0.25) ? 0 : 1;

        // By default the code is set up for running on the local machine. This only changes if run_location == CARC
        bool onLocal = true;
        // For testing purposes using simple saveFolder name
        std::string saveFolder = folder ; // + "/force_check/" + "/met_" + std::to_string(met) +"/" + txLabels[tx] +"/cd8_prolif_" + parameter_levels[prolif_id] + "/cd8_death_" +parameter_levels[death_id] + "/cd8_rec" + parameter_levels[rec_id]; // update to also have the three parameters im sweeping over
        // std::string saveFolder = folder + "/met_" + std::to_string(met) + "/" + txLabels[tx] +"/cd8_prolif_" + parameter_levels[prolif_id] + "/cd8_death_" +parameter_levels[death_id] + "/cd8_rec" + parameter_levels[rec_id]; // update to also have the three parameters im sweeping over
        std::string saveFolderPath = "../../" + saveFolder;
        std::string str  = "python staticParams.py "+ saveFolderPath+" "+replicate_number + " "+ pST + " " + dp_fac + " " + kp_fac; // don't need to run in conda

        if (run_location == "CARC") {
            onLocal = false;
            saveFolderPath = "./" + saveFolder;
            str = "conda run -n bc_env_new python3 genParams.py "+ saveFolderPath+" "+replicate_number + " "+ pST + " " + dp_fac + " " + kp_fac;
        }

        const char *command = str.c_str();
        std::system(command);

        auto start = std::chrono::high_resolution_clock::now();

        Environment model(saveFolderPath, replicate_number,cd8_prolif,cd8_death,cd8_rec,repNum); //can replace with a directory representing any other phenotype state

        model.simulate(1,tx,met,binding_rate_pd1_drug,onLocal);

        auto stop = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> ms_double = stop - start;
        std::cout << "Duration: " << ms_double.count() << std::endl;

    return 0;
}
