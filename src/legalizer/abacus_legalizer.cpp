#include "abacus_legalizer.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>

bool AbacusLegalizer::legalize(PlacementData& data) {
    // Phase 2: Main Legalization Algorithm (Algorithm 1 from paper)
    // Try both ascending and descending order as mentioned in the paper
    
    // Save original state
    PlacementData best_data = data;
    bool found_valid_solution = false;
    double best_total_displacement = std::numeric_limits<double>::infinity();
    
    // Try ascending order first
    for (int sort_dir = 0; sort_dir < 2; ++sort_dir) {
        PlacementData current_data = data;
        bool success = true;
        
        // Step 1: Sort cells by x coordinate
        std::vector<int> cell_order = sort_cells_by_x(current_data.cells, sort_dir == 0);
        
        // Step 2: Legalize each cell one by one
        for (int cell_idx : cell_order) {
            double best_cost = std::numeric_limits<double>::infinity();
            int best_sub_row_idx = -1;
            
            // Get sub-rows sorted by distance from cell
            std::vector<int> sorted_sub_rows = get_sorted_sub_rows(current_data.cells[cell_idx], 
                                                                   current_data.all_sub_rows);
            
            // Try placing cell in each sub-row
            bool found_valid_placement = false;
            
            for (int sub_row_idx : sorted_sub_rows) {
                SubRow* sub_row = current_data.all_sub_rows[sub_row_idx];
                
                // Calculate lower bound cost (vertical movement only)
                double lower_bound = calculate_lower_bound_cost(current_data.cells[cell_idx], *sub_row);
                
                // Early pruning: skip if lower bound exceeds best cost
                if (!should_continue_search(lower_bound, best_cost)) {
                    continue;
                }
                
                // Try to place cell in this sub-row
                double cost;
                bool valid = try_place_cell_in_sub_row(cell_idx, sub_row, current_data, cost);
                
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
                SubRow* best_sub_row = current_data.all_sub_rows[best_sub_row_idx];
                best_sub_row->add_cell(cell_idx, current_data.cells);
                
                // Final placement with PlaceRow
                if (!PlaceRow::place_row(best_sub_row, current_data.cells, 
                                        current_data.max_displacement_constraint, false)) {
                    success = false;
                    break;
                }
                
                current_data.cells[cell_idx].is_legalized = true;
            } else {
                // No valid placement found for this cell
                success = false;
                break;
            }
        }
        
        // Check if this solution is better
        if (success) {
            double total_disp = current_data.calculate_total_displacement();
            if (total_disp < best_total_displacement) {
                best_total_displacement = total_disp;
                best_data = current_data;
                found_valid_solution = true;
            }
        }
    }
    
    if (!found_valid_solution) {
        std::cerr << "Error: No valid placement found" << std::endl;
        return false;
    }
    
    // Apply best solution
    data = best_data;
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
    // Check if cell can fit in sub-row before trying
    // Calculate actual sites needed by considering site alignment
    double total_sites_needed = 0;
    
    // For existing cells in sub-row
    if (!sub_row->cells.empty()) {
        std::vector<int> temp_cells = sub_row->cells;
        std::sort(temp_cells.begin(), temp_cells.end(),
                  [&data](int a, int b) {
                      return data.cells[a].original_x < data.cells[b].original_x;
                  });
        
        double current_x = sub_row->start_x;
        for (int idx : temp_cells) {
            // Each cell starts at a site boundary
            double cell_end = current_x + data.cells[idx].width;
            // Next cell must start at the next site after this cell ends
            double next_site_offset = std::ceil((cell_end - sub_row->start_x) / sub_row->site_width);
            current_x = sub_row->start_x + next_site_offset * sub_row->site_width;
        }
        
        // Add the new cell
        double new_cell_end = current_x + data.cells[cell_index].width;
        total_sites_needed = std::ceil((new_cell_end - sub_row->start_x) / sub_row->site_width);
    } else {
        // Just the new cell
        total_sites_needed = std::ceil(static_cast<double>(data.cells[cell_index].width) / sub_row->site_width);
    }
    
    // Check if total sites needed exceeds available sites
    double available_sites = (sub_row->end_x - sub_row->start_x) / sub_row->site_width;
    if (total_sites_needed > available_sites + PlaceRow::EPSILON) {
        cost = std::numeric_limits<double>::infinity();
        return false;
    }
    
    // Add cell to sub-row
    sub_row->add_cell(cell_index, data.cells);
    
    // Save original positions of all cells in this sub-row
    auto original_positions = save_sub_row_positions(sub_row, data.cells);
    
    // Try to place using PlaceRow (not trial mode - we need actual positions)
    bool valid = PlaceRow::place_row(sub_row, data.cells, 
                                    data.max_displacement_constraint, false);
    
    if (valid) {
        // Verify no overlaps after placement
        std::vector<int> cell_indices = sub_row->cells;
        std::sort(cell_indices.begin(), cell_indices.end(),
                  [&data](int a, int b) {
                      return data.cells[a].current_x < data.cells[b].current_x;
                  });
        
        // Check for overlaps
        for (int i = 1; i < static_cast<int>(cell_indices.size()); ++i) {
            double prev_end = data.cells[cell_indices[i-1]].current_x + 
                             data.cells[cell_indices[i-1]].width;
            double curr_start = data.cells[cell_indices[i]].current_x;
            
            if (prev_end > curr_start + PlaceRow::EPSILON) {
                valid = false;
                break;
            }
        }
        
        // Check boundaries
        if (valid) {
            for (int idx : cell_indices) {
                if (data.cells[idx].current_x < sub_row->start_x - PlaceRow::EPSILON ||
                    data.cells[idx].current_x + data.cells[idx].width > 
                    sub_row->end_x + PlaceRow::EPSILON) {
                    valid = false;
                    break;
                }
            }
        }
    }
    
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