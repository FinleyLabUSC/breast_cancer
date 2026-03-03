#include <array>

#include "../inc/Environment.h"
#include "../inc/CellGrid.h"

// A data structure for quick access to cells in the neighborhood of other cells
CellGrid::CellGrid(double width, double max_dist)
{
    grid_width = width; // Width of bins
    max_d = max_dist; // Max distance for neighbors / influences

    // We assume there won't be more than 256 grid cells to begin (16*16)
    cell_grid.reserve(256);
    cancer_locs.reserve(256);
}

// Given a list of shared pointers to cells, update the state of the grid
void CellGrid::update_grid(std::vector<std::shared_ptr<RS_Cell>>& cell_list)
{
    // update_grid() assumes that the grid has ALREADY been cleared
    for (int i = 0; i < cell_list.size(); ++i)
    {
        // Get bin location & bin key (from Szudzik elegant pairing function)
        auto loc = cell_list[i]->x;
        int x_bin = std::floor(loc[0] / grid_width);
        int y_bin = std::floor(loc[1] / grid_width);
        int key = get_szudzik(x_bin, y_bin);

        // If x is cancer, add it to cancer_locs; this will prevent multi-inserts
        cancer_locs.insert(key);

        auto new_shared_ptr(cell_list[i]);
        if (cell_grid.find(key) == cell_grid.end())
        {   // The bin does not exist
            std::vector<std::shared_ptr<RS_Cell>> temp_vec;
            temp_vec.reserve(512); // Generous limit on # cells per grid-square based on NK radius
            temp_vec.push_back(new_shared_ptr);
            cell_grid[key] = temp_vec;
            continue; // The shared_ptr has been allocated, so we can skip
        }

        // The bin exists
        cell_grid[key].push_back(new_shared_ptr);
    }
}

// Get the nearest neighbor of ANY cell type, given the grid's max distance
std::array<double, 2> CellGrid::get_NN(std::array<double, 2> loc)
{
    int x_bin = std::floor(loc[0] / grid_width);
    int y_bin = std::floor(loc[1] / grid_width);

    bool found_neighbor = false;
    double curr_nn_dsq = 100000; // squared distance to current NN
    std::array<double, 2> nn_loc = {0,0};
    int max_ring = std::ceil(max_d/grid_width);

    // This loops searches increasing rings up to the max_ring
    for (int r = 1; r <= max_ring; ++r)
    {   // This loop will search the (2r+1)x(2r+1) grid around the central box
        for (int i = -r; i <= r; ++i)
        {
            for (int j = -r; j <= r; ++j)
            {
                if (r != 1 && (j != r) && (i != r)){continue;} // If r != 1, then we only want to search the outer ring
                int key = get_szudzik(x_bin + i, y_bin + j); // Generate the access key
                if (cell_grid.find(key) == cell_grid.end()){continue;} // Skip if the bin doesn't exist
                for (auto cell : cell_grid[key])
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

// Get the nearest neighbor of a SPECIFIC cell type, within a certain max distance
std::array<double, 2> CellGrid::get_filter_NN(std::array<double, 2> loc, double max_dist, int state)
{
    // Location for neighbors search
    int x_bin = std::floor(loc[0] / grid_width);
    int y_bin = std::floor(loc[1] / grid_width);

    bool found_neighbor = false;
    double curr_nn_dsq = 100000; // squared distance to current NN
    std::array<double, 2> nn_loc = {0,0};
    int max_ring = std::ceil(max_dist/grid_width); // use max_dist instead of max_d

    // This loops searches increasing rings up to the max_ring
    for (int r = 1; r <= max_ring; ++r)
    {   // This loop will search the (2r+1)x(2r+1) grid around the central box
        for (int i = -r; i <= r; ++i)
        {
            for (int j = -r; j <= r; ++j)
            {
                if (r != 1 && (j != r) && (i != r)){continue;} // If r != 1, then we only want to search the outer ring
                int key = get_szudzik(x_bin + i, y_bin + j); // Generate the access key
                if (cell_grid.find(key) == cell_grid.end()){continue;} // Skip if the bin doesn't exist
                for (auto cell : cell_grid[key])
                {
                    if (cell->state != state) {continue;} // skip if cell isn't in the right state

                    // check distance between cells
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
        if (found_neighbor == true && curr_nn_dsq <= r*grid_width*r*grid_width && curr_nn_dsq <= max_dist*max_dist)
        {
            return nn_loc;
        }
    }
    // In the case where there is no NN within the correct distance, return the cell's location
    return loc;
}

// Gets all cell neighbors within a max_d specified by the grid (usually 100 microns)
std::vector<int> CellGrid::get_neighbors(std::array<double, 2> loc, int runtime_idx)
{
    std::vector<int> neighbors; // place to store neighbors
    int x_bin = std::floor(loc[0] / grid_width); // x_bin
    int y_bin = std::floor(loc[1] / grid_width); // y_bin

    // Iterate over the 3x3 grid around this bin; check to make sure bins exist
    for (int i = x_bin - 1; i <= x_bin + 1; ++i)
    {
        for (int j = y_bin - 1; j <= y_bin + 1; ++j)
        {
            int key = get_szudzik(i, j);
            if (cell_grid.find(key) == cell_grid.end()){continue;}
            for (auto cell : cell_grid[key])
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

int CellGrid::get_szudzik(int x, int y)
{
    // Create unique int from two ints:
    int new_x = (x < 0) ? -2*x - 1 : 2*x; // Map x_bin to +Z
    int new_y = (y < 0) ? -2*y - 1 : 2*y; // Map y_bin to +Z
    return (new_x >= new_y) ? new_x*new_x+new_x+new_y : new_y*new_y+new_x; // Szudzik elegant pairing function
}

void CellGrid::clear_grid()
{
    cell_grid.clear();
}

// Clear and reassemble the cell_grid and cancer_grid
void Environment::update_grids()
{
    cell_grid.clear_grid();
    cell_grid.update_grid(cell_list);
}
