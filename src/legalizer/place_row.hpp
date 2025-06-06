// Simplified place_row.hpp - matching reference implementation style
#ifndef PLACE_ROW_HPP
#define PLACE_ROW_HPP

#include "../data_structure/data_structure.hpp"
#include <vector>
#include <utility>

class PlaceRow {
public:
    // Trial mode - returns (isValid, cost) pair like reference
    // isValid: -1 if failed, >= 0 if successful
    // cost: displacement of the new cell
    static std::pair<int, double> place_row_trial(SubRow* sub_row, 
                                                  Cell* new_cell,
                                                  const std::vector<Cell>& cells,
                                                  double max_displacement_constraint,
                                                  bool add_penalty);
    
    // Final mode - actually modifies cluster state
    static void place_row_final(SubRow* sub_row, 
                               Cell* new_cell,
                               std::vector<Cell>& cells);
    
    // Get site-aligned x position
    static double get_site_x(double x, double min_x, int site_width);
};

#endif // PLACE_ROW_HPP