// abacus_legalizer.hpp - Modified version
#ifndef ABACUS_LEGALIZER_HPP
#define ABACUS_LEGALIZER_HPP

#include "../data_structure/data_structure.hpp"
#include "../legalizer/place_row.hpp"
#include <vector>
#include <limits>
#include <algorithm>

class AbacusLegalizer {
public:
    // Main legalization function
    bool legalize(PlacementData& data);

private:
    // Sort cells by their original x position
    std::vector<int> sort_cells_by_x(const std::vector<Cell>& cells, bool ascending = true);
    
    // Find the nearest row to a cell
    int find_nearest_row(const Cell& cell, const PlacementData& data);
    
    // Calculate lower bound cost for placing a cell in a sub-row (vertical movement only)
    double calculate_lower_bound_cost(const Cell& cell, const SubRow& sub_row);
    
    // Try to place a cell in a specific sub-row using trial mode
    PlacementResult try_place_cell_trial(int cell_index, SubRow* sub_row, 
                                       const PlacementData& data, 
                                       double current_best_cost,
                                       bool add_penalty = true);  // Added add_penalty parameter
    
    // Apply site alignment to all cells after placement
    void apply_site_alignment_all(PlacementData& data);
    
    // Get the index of a sub-row in the all_sub_rows vector
    int get_sub_row_index(SubRow* sub_row, const PlacementData& data);
    
    // Get sub-rows for a specific row
    std::vector<SubRow*> get_sub_rows_for_row(const PlacementData& data, int row_idx);
    
    // Find best sub-row within a row
    std::pair<SubRow*, PlacementResult> find_best_sub_row_in_row(
        int cell_idx, int row_idx, const PlacementData& data, 
        double current_best_cost, bool add_penalty);
};

#endif // ABACUS_LEGALIZER_HPP