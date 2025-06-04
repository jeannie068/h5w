#ifndef ABACUS_LEGALIZER_HPP
#define ABACUS_LEGALIZER_HPP

#include "../data_structure/data_structure.hpp"
#include "place_row.hpp" // Ensure this is the new place_row.hpp
#include <vector>
#include <limits>
#include <algorithm>

class AbacusLegalizer {
public:
    bool legalize(PlacementData& data);

private:
    std::vector<int> sort_cells_by_x(const std::vector<Cell>& cells, bool ascending = true);
    int find_nearest_row_idx(const Cell& cell, const PlacementData& data);
    
    // Tries to place cell_index into sub_row using PlaceRow::place_row_trial
    PlacementTrialResult try_place_cell_in_sub_row_trial(
        const SubRow* sub_row,
        int cell_idx,
        const PlacementData& placement_data,
        bool check_penalty_for_all_cluster_cells);

    // Finds the best sub-row within a given row for a cell
    std::pair<SubRow*, PlacementTrialResult> find_best_sub_row_in_row(
        int cell_idx, 
        int row_idx, 
        PlacementData& placement_data, // Not const, as subrow might be modified in trial if not careful (but shouldn't)
                                       // Making it const to enforce no modification in trial search
        const PlacementData& const_placement_data, // For trial calls
        double current_best_overall_cost,
        bool check_penalty_for_all_cluster_cells);
    
    // Applies site alignment to all cells in all sub-rows.
    // This might be redundant if place_row_final's determine_and_apply_final_positions handles it.
    // Kept for now as a final sweep if needed.
    void apply_site_alignment_to_all_subrows(PlacementData& data);
};

#endif // ABACUS_LEGALIZER_HPP
