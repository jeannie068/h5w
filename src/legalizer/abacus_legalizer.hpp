// Simplified abacus_legalizer.hpp
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
    
    // Try to place a cell in a specific row, returns (sub_row_idx, cost)
    std::pair<int, double> try_place_in_row(Cell* cell, int row_idx, 
                                            const PlacementData& data, 
                                            bool add_penalty);
    
    // Determine final positions for all cells (site alignment)
    void determine_final_positions(PlacementData& data);
};

#endif // ABACUS_LEGALIZER_HPP