#include "abacus_legalizer.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>

bool AbacusLegalizer::legalize(PlacementData& data) {
    // Phase 2: Main Legalization Algorithm (Algorithm 1 from paper)
    
    // Step 1: Sort cells by x coordinate (ascending order)
    std::vector<int> cell_order = sort_cells_by_x(data.cells, true);
    
    // Step 2: Legalize each cell one by one
    for (int cell_idx : cell_order) {
        double best_cost = std::numeric_limits<double>::infinity();
        int best_sub_row_idx = -1;
        
        // Get sub-rows sorted by distance from cell
        std::vector<int> sorted_sub_rows = get_sorted_sub_rows(data.cells[cell_idx], 
                                                               data.all_sub_rows);
        
        // Try placing cell in each sub-row
        bool found_valid_placement = false;
        
        for (int sub_row_idx : sorted_sub_rows) {
            SubRow* sub_row = data.all_sub_rows[sub_row_idx];
            
            // Calculate lower bound cost (vertical movement only)
            double lower_bound = calculate_lower_bound_cost(data.cells[cell_idx], *sub_row);
            
            // Early pruning: skip if lower bound exceeds best cost
            if (!should_continue_search(lower_bound, best_cost)) {
                // If searching upward/downward and lower bound exceeds best cost,
                // we can stop searching in that direction
                continue;
            }
            
            // Try to place cell in this sub-row
            double cost;
            bool valid = try_place_cell_in_sub_row(cell_idx, sub_row, data, cost);
            
            // Remove cell after trial (it was added in try_place_cell_in_sub_row)
            sub_row->remove_cell(cell_idx);
            
            if (valid && cost < best_cost) {
                best_cost = cost;
                best_sub_row_idx = sub_row_idx;
                found_valid_placement = true;
            }
        }
        
        // Place cell in best sub-row (if found)
        if (found_valid_placement && best_sub_row_idx >= 0) {
            SubRow* best_sub_row = data.all_sub_rows[best_sub_row_idx];
            best_sub_row->add_cell(cell_idx, data.cells);
            
            // Final placement with PlaceRow
            if (!PlaceRow::place_row(best_sub_row, data.cells, 
                                    data.max_displacement_constraint, false)) {
                std::cerr << "Error: Failed to place cell " << data.cells[cell_idx].name 
                         << " in best sub-row" << std::endl;
                return false;
            }
            
            data.cells[cell_idx].is_legalized = true;
        } else {
            // No valid placement found for this cell
            std::cerr << "Error: No valid placement found for cell " 
                     << data.cells[cell_idx].name << std::endl;
            return false;
        }
    }
    
    // Verify all cells are legalized
    for (const auto& cell : data.cells) {
        if (!cell.is_legalized) {
            std::cerr << "Error: Cell " << cell.name << " is not legalized" << std::endl;
            return false;
        }
    }
    
    return true;
}

std::vector<int> AbacusLegalizer::sort_cells_by_x(const std::vector<Cell>& cells, 
                                                  bool ascending) {
    std::vector<int> indices(cells.size());
    for (int i = 0; i < static_cast<int>(cells.size()); ++i) {
        indices[i] = i;
    }
    
    if (ascending) {
        std::sort(indices.begin(), indices.end(),
                  [&cells](int a, int b) {
                      return cells[a].original_x < cells[b].original_x;
                  });
    } else {
        std::sort(indices.begin(), indices.end(),
                  [&cells](int a, int b) {
                      return cells[a].original_x > cells[b].original_x;
                  });
    }
    
    return indices;
}

double AbacusLegalizer::calculate_lower_bound_cost(const Cell& cell, 
                                                  const SubRow& sub_row) {
    // Lower bound = vertical movement only
    double dy = sub_row.y - cell.original_y;
    return std::abs(dy);
}

double AbacusLegalizer::calculate_actual_cost(const Cell& cell, 
                                            double new_x, double new_y) {
    double dx = new_x - cell.original_x;
    double dy = new_y - cell.original_y;
    return std::sqrt(dx * dx + dy * dy);
}

bool AbacusLegalizer::try_place_cell_in_sub_row(int cell_index, SubRow* sub_row,
                                               PlacementData& data, double& cost) {
    // Add cell to sub-row
    sub_row->add_cell(cell_index, data.cells);
    
    // Save original positions of all cells in this sub-row
    auto original_positions = save_sub_row_positions(sub_row, data.cells);
    
    // Try to place using PlaceRow (not trial mode - we need actual positions)
    bool valid = PlaceRow::place_row(sub_row, data.cells, 
                                    data.max_displacement_constraint, false);
    
    if (valid) {
        // Calculate cost for this placement using actual position from PlaceRow
        cost = calculate_actual_cost(data.cells[cell_index], 
                                   data.cells[cell_index].current_x, 
                                   data.cells[cell_index].current_y);
    } else {
        cost = std::numeric_limits<double>::infinity();
    }
    
    // Restore original positions (since this is just a trial)
    restore_sub_row_cells(sub_row, data.cells, original_positions);
    
    return valid;
}

std::vector<int> AbacusLegalizer::get_sorted_sub_rows(const Cell& cell,
                                                     const std::vector<SubRow*>& all_sub_rows) {
    std::vector<std::pair<double, int>> distances;
    
    for (int i = 0; i < static_cast<int>(all_sub_rows.size()); ++i) {
        // Calculate distance to sub-row (use center of sub-row)
        double sub_row_center_x = (all_sub_rows[i]->start_x + all_sub_rows[i]->end_x) / 2.0;
        double dx = sub_row_center_x - cell.original_x;
        double dy = all_sub_rows[i]->y - cell.original_y;
        double dist = std::sqrt(dx * dx + dy * dy);
        distances.push_back({dist, i});
    }
    
    // Sort by distance
    std::sort(distances.begin(), distances.end());
    
    std::vector<int> sorted_indices;
    for (const auto& pair : distances) {
        sorted_indices.push_back(pair.second);
    }
    
    return sorted_indices;
}

bool AbacusLegalizer::should_continue_search(double lower_bound, double best_cost) {
    return lower_bound < best_cost;
}

void AbacusLegalizer::restore_sub_row_cells(SubRow* sub_row, std::vector<Cell>& cells,
                                          const std::vector<std::pair<double, double>>& original_positions) {
    for (int i = 0; i < static_cast<int>(sub_row->cells.size()) && 
                    i < static_cast<int>(original_positions.size()); ++i) {
        int cell_idx = sub_row->cells[i];
        cells[cell_idx].current_x = original_positions[i].first;
        cells[cell_idx].current_y = original_positions[i].second;
    }
}

std::vector<std::pair<double, double>> AbacusLegalizer::save_sub_row_positions(
    const SubRow* sub_row, const std::vector<Cell>& cells) {
    std::vector<std::pair<double, double>> positions;
    
    for (int cell_idx : sub_row->cells) {
        positions.push_back({cells[cell_idx].current_x, cells[cell_idx].current_y});
    }
    
    return positions;
}