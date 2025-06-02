#include "../inc/ModelUtil.h"

// TODO this gets removed anyway, so I've not worrying about it right now.
int getRandomNumber(int maxNumber) {

    //takes in an integer maxNumber and returns a random integer from 0 to maxNumber-1
    std::mt19937 gen(1); // Mersenne Twister 32-bit PRNG using seed from random device
    std::uniform_int_distribution<> distrib(0, maxNumber-1); // Define the range

    // return distrib(gen); // Generate and return a random number within the specified range
    return 0;
}

/* Takes in a 2D vector of strings (@param vec), and extracts a specific row from it, specified by @param row_idx*/
std::vector<std::string> get2dvecrow(std::vector<std::vector<std::string> >& vec, size_t row_idx){
    std::vector<std::string> res; 
    for(size_t i = 0; i < vec[row_idx].size(); ++i){
        
        res.push_back(vec[row_idx][i]); 
    }

    return res; 

}