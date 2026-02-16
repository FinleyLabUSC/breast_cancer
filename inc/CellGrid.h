#ifndef BREAST_CANCER_CELLGRID_H
#define BREAST_CANCER_CELLGRID_H

#include "RS_Cell.h"

class  CellGrid
{
public:
    // virtual destructor
    virtual ~CellGrid() = default;

    // members
    CellGrid(double width, double max_dist);
    std::unordered_map<int, std::unordered_map<int, std::vector<std::shared_ptr<RS_Cell>>>> cell_grid;
    double grid_width = 100;
    double max_d = 100;

    // fuctions
    virtual void update_grid(std::vector<std::shared_ptr<RS_Cell>>& cell_list);
    void clear_grid();
    std::array<double, 2> get_NN(std::array<double, 2> loc);
    std::vector<int> get_neighbors(std::array<double, 2> loc, int runtime_idx);
};

class Restricted_CellGrid final : public CellGrid
{
public:
    // slightly different initialization
    Restricted_CellGrid(double width, double max_dist, int type);
    int filter_type = -1;

    // only change the update_grid function
    void update_grid(std::vector<std::shared_ptr<RS_Cell>>& cell_list) override;
};

#endif //BREAST_CANCER_CELLGRID_H