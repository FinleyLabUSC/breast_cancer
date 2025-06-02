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
    std::string folder = argv[1];
    std::string set = argv[2];
    std::string pST = argv[3];
    std::string dp_fac = argv[4];
    std::string kp_fac = argv[5];
    int tx = std::stoi(argv[6]);
    int met = std::stoi(argv[7]);
    int number_of_chemo_doses = std::stoi(argv[8]);
    int repNum = std::stoi(argv[2]);

    std::vector<std::string> txLabels = {"control","chemo","combo"};
    std::vector<std::string> metLabels = {"met1","met2"};
    // TODO you probably need to change this for the cluster
    std::string saveFolder = "../../" + folder + "/" + metLabels[met-1] + "/" + txLabels[tx] + "/dose_" + std::to_string(number_of_chemo_doses);

    //std::string str = "rm -r ./"+folder+"/set_" + set;
    //const char *command = str.c_str();
    //std::system(command);

    // str = "python genParams.py ./"+folder+"/set_"+set+" "+set;
    std::string str = "conda run -n bc_env python3 ../genParams.py "+ saveFolder+" "+set + " "+ pST + " " + dp_fac + " " + kp_fac;
    const char *command = str.c_str();
    std::system(command);

    auto start = std::chrono::high_resolution_clock::now();

    Environment model(saveFolder, set, "../phenotype_out/",repNum); //can replace with a directory representing any other phenotype state

   // model.generateNums();
    model.simulate(1,tx,met,number_of_chemo_doses);

    auto stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> ms_double = stop - start;
    std::cout << "Duration: " << ms_double.count() << std::endl;
    return 0;
}
