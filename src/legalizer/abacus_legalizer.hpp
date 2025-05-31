#ifndef ABACUS_LEGALIZER_HPP
#define ABACUS_LEGALIZER_HPP

#include "../data_structure/data_structure.hpp"
#include "../legalizer/place_row.hpp"
#include <vector>
#include <limits>

class AbacusLegalizer {
public:
    // Main legalization function
    bool legalize(PlacementData& data);

private:
    // Sort cells by their original x position
    std::vector<int> sort_cells_by_x(const std::vector<Cell>& cells, bool ascending = true);
    
    // Calculate lower bound cost for placing a cell in a sub-row (vertical movement only)
    double calculate_lower_bound_cost(const Cell& cell, const SubRow& sub_row);
    
    // Calculate actual cost (Euclidean displacement) for a cell
    double calculate_actual_cost(const Cell& cell, double new_x, double new_y);
    
    // Try to place a cell in a specific sub-row
    // Returns true if placement is valid, false otherwise
    bool try_place_cell_in_sub_row(int cell_index, SubRow* sub_row, 
                                   PlacementData& data, double& cost);
    
    // Get sub-rows sorted by distance from cell's original position
    std::vector<int> get_sorted_sub_rows(const Cell& cell, 
                                         const std::vector<SubRow*>& all_sub_rows);
    
    // Check if we should continue searching rows in a given direction
    bool should_continue_search(double lower_bound, double best_cost);
    
    // Restore original cell positions for cells in a sub-row
    void restore_sub_row_cells(SubRow* sub_row, std::vector<Cell>& cells,
                              const std::vector< std::pair<double, double> >& original_positions);
    
    // Save current positions of cells in a sub-row
    std::vector< std::pair<double, double> > save_sub_row_positions(const SubRow* sub_row, 
                                                                  const std::vector<Cell>& cells);
};

#endif // ABACUS_LEGALIZER_HPP