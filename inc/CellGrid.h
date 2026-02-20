#ifndef BREAST_CANCER_CELLGRID_H
#define BREAST_CANCER_CELLGRID_H

#include "RS_Cell.h"
#include <unordered_set>

class  CellGrid
{
public:
    // virtual destructor
    virtual ~CellGrid() = default;

    // members
    CellGrid(double width, double max_dist);
    std::unordered_map<int, std::vector<std::shared_ptr<RS_Cell>>> cell_grid;
    std::unordered_set<int> cancer_locs;
    double grid_width = 100;
    double max_d = 100;

    // fuctions
    virtual void update_grid(std::vector<std::shared_ptr<RS_Cell>>& cell_list);
    void clear_grid();
    std::array<double, 2> get_NN(std::array<double, 2> loc);
    std::array<double, 2> get_filter_NN(std::array<double, 2> loc, double max_dist, int state);
    std::vector<int> get_neighbors(std::array<double, 2> loc, int runtime_idx);
    int get_szudzik(int x, int y);
};


#endif //BREAST_CANCER_CELLGRID_H