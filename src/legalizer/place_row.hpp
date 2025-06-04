#ifndef PLACE_ROW_HPP
#define PLACE_ROW_HPP

#include "../data_structure/data_structure.hpp"
#include <vector>
#include <limits> // Required for std::numeric_limits
#include <memory> // Required for std::shared_ptr

// Structure to hold placement result for trial mode
struct PlacementTrialResult {
    bool valid;
    double cost;                          // Displacement of the new cell
    double optimal_x_for_new_cell;      // Calculated optimal X for the new cell being trialed

    PlacementTrialResult() : valid(false), cost(std::numeric_limits<double>::infinity()), optimal_x_for_new_cell(0.0) {}
};

class PlaceRow {
public:
    static constexpr double EPSILON = 1e-9;

    // Trial mode: Calculates the cost and optimal X for the new_cell_index if it were placed in sub_row.
    // Does NOT modify sub_row or global cell data.
    static PlacementTrialResult place_row_trial(
        const SubRow* sub_row, 
        int new_cell_index,
        const PlacementData& placement_data,
        double max_displacement_constraint,
        bool check_penalty_for_all_cluster_cells);

    // Final mode: Actually places the new_cell_index into the sub_row.
    // Modifies sub_row->last_cluster and updates positions of affected cells in placement_data.
    // Also adds new_cell_index to sub_row->cells.
    static void place_row_final(
        SubRow* sub_row,
        int new_cell_index,
        PlacementData& placement_data);
    
private:
    // Helper to determine and apply final positions for all cells in a sub-row based on its cluster chain.
    // Called by place_row_final. Handles site alignment.
    static void determine_and_apply_final_positions(SubRow* sub_row, PlacementData& placement_data);
};

#endif // PLACE_ROW_HPP
