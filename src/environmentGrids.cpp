#include <array>

#include "../inc/Environment.h"
#include "../inc/CellGrid.h"

CellGrid::CellGrid(double width, double max_dist)
{
    grid_width = width; // Width of bins
    max_d = max_dist; // Max distance for neighbors / influences
}

void CellGrid::update_grid(std::vector<std::shared_ptr<RS_Cell>>& cell_list)
{
    // update_grid() assumes that the grid has ALREADY been cleared
    for (int i = 0; i < cell_list.size();++i)
    {
        // Get bin location
        auto loc = cell_list[i]->x;
        int x_bin = std::floor(loc[0] / grid_width);
        int y_bin = std::floor(loc[1] / grid_width);

        // First see if the x-row exists
        if (cell_grid.find(x_bin) == cell_grid.end())
        {
            // If the x-row doesn't exist, make it & allocate the y-bin to it (no y-bins exist yet)
            std::unordered_map<int, std::vector<std::shared_ptr<RS_Cell>>> temp_map;
            auto new_shared_ptr(cell_list[i]);
            std::vector temp_vect = {new_shared_ptr};
            temp_map[y_bin] = temp_vect;
            cell_grid[x_bin] = temp_map;
            continue; // The shared_ptr has been allocated, so we can skip
        }

        // The x-row exists, check to see if the y-bin exists
        std::unordered_map<int, std::vector<std::shared_ptr<RS_Cell>>> x_row = cell_grid[x_bin];
        if (x_row.find(y_bin) == x_row.end())
        {
            // If the y-bin doesn't exist, make it & allocate the shared_ptr to it
            auto new_shared_ptr(cell_list[i]);
            std::vector temp_vect = {new_shared_ptr};
            cell_grid[x_bin][y_bin] = temp_vect;
            continue; // The shared_ptr has been allocated, so we can skip
        }

        // The x-row AND the y-bin exists, so we just need to extend the vector
        auto new_shared_ptr(cell_list[i]);
        cell_grid[x_bin][y_bin].push_back(new_shared_ptr);
    }
}

void CellGrid::clear_grid()
{
    cell_grid.clear();
}

Restricted_CellGrid::Restricted_CellGrid(double width, double max_dist, int type) :
    CellGrid(width, max_dist)
{
    filter_type = type; // Celltype to keep
    grid_width = width; // Width of bins
}

void Restricted_CellGrid::update_grid(std::vector<std::shared_ptr<RS_Cell>>& cell_list)
{
    // update_grid() assumes that the grid has ALREADY been cleared
    for (int i = 0; i < cell_list.size();++i)
    {
        // Check cell STATE status
        if (cell_list[i]->state != filter_type){continue;}

        // Get bin location
        auto loc = cell_list[i]->x;
        int x_bin = std::floor(loc[0] / grid_width);
        int y_bin = std::floor(loc[1] / grid_width);

        // First see if the x-row exists
        if (cell_grid.find(x_bin) == cell_grid.end())
        {
            // If the x-row doesn't exist, make it & allocate the y-bin to it (no y-bins exist yet)
            std::unordered_map<int, std::vector<std::shared_ptr<RS_Cell>>> temp_map;
            auto new_shared_ptr(cell_list[i]);
            std::vector temp_vect = {new_shared_ptr};
            temp_map[y_bin] = temp_vect;
            cell_grid[x_bin] = temp_map;
            continue; // The shared_ptr has been allocated, so we can skip
        }

        // The x-row exists, check to see if the y-bin exists
        std::unordered_map<int, std::vector<std::shared_ptr<RS_Cell>>> x_row = cell_grid[x_bin];
        if (x_row.find(y_bin) == x_row.end())
        {
            // If the y-bin doesn't exist, make it & allocate the shared_ptr to it
            auto new_shared_ptr(cell_list[i]);
            std::vector temp_vect = {new_shared_ptr};
            cell_grid[x_bin][y_bin] = temp_vect;
            continue; // The shared_ptr has been allocated, so we can skip
        }

        // The x-row AND the y-bin exists, so we just need to extend the vector
        auto new_shared_ptr(cell_list[i]);
        cell_grid[x_bin][y_bin].push_back(new_shared_ptr);
    }
}

std::vector<int> CellGrid::get_neighbors(std::array<double, 2> loc, int runtime_idx)
{
    std::vector<int> neighbors; // place to store neighbors
    int x_bin = std::floor(loc[0] / grid_width); // x_bin
    int y_bin = std::floor(loc[1] / grid_width); // y_bin

    // Iterate over the 3x3 grid around this bin; check to make sure bins exist
    for (int i = x_bin - 1; i <= x_bin + 1; ++i)
    {
        if (cell_grid.find(i) == cell_grid.end()){continue;}
        for (int j = y_bin - 1; j <= y_bin + 1; ++j)
        {
            if (cell_grid[i].find(j) == cell_grid[i].end()){continue;}
            for (auto cell : cell_grid[i][j])
            {   // check distance between cells
                std::array<double, 2> this_loc = cell->x;
                double dx = this_loc[0] - loc[0];
                double dy = this_loc[1] - loc[1];
                if (dx*dx + dy*dy <= max_d*max_d && cell->runtime_index != runtime_idx)
                {   // push back the current runtime_index if it's not the cell's
                    neighbors.push_back(cell->runtime_index);
                }
            }
        }
    }

    return neighbors;
}

std::array<double, 2> CellGrid::get_NN(std::array<double, 2> loc)
{
    int x_bin = std::floor(loc[0] / grid_width);
    int y_bin = std::floor(loc[1] / grid_width);

    bool found_neighbor = false;
    double curr_nn_dsq = 100000; // squared distance to current NN
    std::array<double, 2> nn_loc = {0,0};
    int max_ring = std::round(max_d/grid_width);

    for (int r = 1; r <= max_ring; ++r)
    {   // This loop will search the (2r+1)x(2r+1) grid around the central box
        for (int i = -r; i <= r; ++i)
        {
            if (cell_grid.find(x_bin + i) == cell_grid.end()){continue;}
            for (int j = -r; j <= r; ++j)
            {
                // If r != 1, then we only want to search the outer ring
                if (r != 1 && (j != r) && (i != r)){continue;}
                // If r == 1 or (j or i == r) (e.g., we are on the outer ring), proceed
                if (cell_grid[x_bin + i].find(y_bin + j) == cell_grid[x_bin + i].end()){continue;}
                for (auto cell : cell_grid[x_bin + i][y_bin + j])
                {   // check distance between cells
                    std::array<double, 2> this_loc = cell->x;
                    double dx = this_loc[0] - loc[0];
                    double dy = this_loc[1] - loc[1];
                    if (dx*dx + dy*dy <= curr_nn_dsq)
                    {
                        nn_loc = cell->x; // store current location
                        curr_nn_dsq = dx*dx + dy*dy; // store squared distance
                        found_neighbor = true; // set flag to true
                    }
                }
            }
        }

        // Every time a loop terminates, we need to check if the NN has been found & is closer than max_d
        // BUT, the nn must be within r*grid_width radius, or we need to search the next ring
        if (found_neighbor == true && curr_nn_dsq <= r*grid_width*r*grid_width && curr_nn_dsq <= max_d*max_d)
        {
            return nn_loc;
        }
    }
    // In the case where there is no NN within the correct distance, return the cell's location
    return loc;
}

// Clear and reassemble the cell_grid and cancer_grid
void Environment::update_grids(bool do_cancer)
{
    cell_grid.clear_grid();
    cell_grid.update_grid(cell_list);
    cancer_grid.clear_grid();
    cancer_grid.update_grid(cell_list);
}
