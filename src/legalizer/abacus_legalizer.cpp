// abacus_legalizer.cpp - Modified to use incremental PlaceRow
#include "abacus_legalizer.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>

// Modified try_place_cell_trial to use incremental placement
PlacementResult AbacusLegalizer::try_place_cell_trial(int cell_index, SubRow* sub_row,
                                                     const PlacementData& data, 
                                                     double current_best_cost,
                                                     bool add_penalty) {
    PlacementResult result;
    result.valid = false;
    result.cost = std::numeric_limits<double>::infinity();
    
    // Check if cell width can fit in sub-row considering site alignment
    int cell_sites_needed = static_cast<int>(std::ceil(static_cast<double>(data.cells[cell_index].width) / 
                                                       sub_row->site_width));
    
    // Calculate current used sites
    int used_sites = 0;
    for (int idx : sub_row->cells) {
        int sites = static_cast<int>(std::ceil(static_cast<double>(data.cells[idx].width) / 
                                               sub_row->site_width));
        used_sites += sites;
    }
    
    int available_sites = static_cast<int>((sub_row->end_x - sub_row->start_x) / sub_row->site_width);
    int remaining_sites = available_sites - used_sites;
    
    // Check if cell can fit in terms of sites
    if (cell_sites_needed > remaining_sites) {
        return result;
    }
    
    // Use incremental trial mode of PlaceRow
    result = PlaceRow::place_row_incremental_trial(sub_row, cell_index, data.cells, 
                                                   data.max_displacement_constraint, add_penalty);
    
    return result;
}

bool AbacusLegalizer::legalize(PlacementData& data) {
    std::cout << "\n===== Starting Abacus Legalization (Incremental) =====" << std::endl;
    
    // Try both ascending and descending order
    PlacementData best_data = data;
    bool found_valid_solution = false;
    double best_total_displacement = std::numeric_limits<double>::infinity();

    int max_dir_num = (data.cells.size() < 220000) ? 2 : 1; // Use both directions for smaller datasets, only one for larger
    
    for (int sort_dir = 0; sort_dir < max_dir_num; ++sort_dir) {
        std::cout << "\n--- Trying sort direction: " << (sort_dir == 0 ? "ascending" : "descending") << " ---" << std::endl;
        
        PlacementData current_data = data;
        
        // Initialize all sub-rows' last_cluster to nullptr
        for (auto& row : current_data.rows) {
            for (auto& sub_row : row.sub_rows) {
                sub_row.last_cluster = nullptr;
                sub_row.cells.clear();
            }
        }
        
        bool success = true;
        
        // Sort cells by x coordinate
        std::vector<int> cell_order = sort_cells_by_x(current_data.cells, sort_dir);
        
        // Legalize each cell one by one
        for (int i = 0; i < static_cast<int>(cell_order.size()); ++i) {
            int cell_idx = cell_order[i];
            
            double best_cost = std::numeric_limits<double>::infinity();
            SubRow* best_sub_row = nullptr;
            PlacementResult best_result;
            
            // Find nearest row based on Y coordinate
            int base_row_idx = find_nearest_row(current_data.cells[cell_idx], current_data);
            bool found_valid_placement = false;
            
            // Try with penalty first, then without if needed
            for (int add_penalty = 1; add_penalty >= 0; --add_penalty) {
                // Search upward from base row
                for (int row_idx = base_row_idx; row_idx >= 0; --row_idx) {
                    // Calculate lower bound (vertical distance only)
                    double vertical_distance = std::abs(current_data.cells[cell_idx].original_y - 
                                                      current_data.rows[row_idx].y);
                    
                    // Pruning: if vertical distance alone exceeds best cost, stop searching upward
                    if (vertical_distance >= best_cost) {
                        break;
                    }
                    
                    // Try sub-rows in this row
                    auto [sub_row, result] = find_best_sub_row_in_row(
                        cell_idx, row_idx, current_data, best_cost, add_penalty);
                    
                    if (result.valid && result.cost < best_cost) {
                        best_cost = result.cost;
                        best_sub_row = sub_row;
                        best_result = result;
                        found_valid_placement = true;
                    }
                }
                
                // Search downward from base row
                for (int row_idx = base_row_idx + 1; row_idx < static_cast<int>(current_data.rows.size()); ++row_idx) {
                    // Calculate lower bound (vertical distance only)
                    double vertical_distance = std::abs(current_data.cells[cell_idx].original_y - 
                                                      current_data.rows[row_idx].y);
                    
                    // Pruning: if vertical distance alone exceeds best cost, stop searching downward
                    if (vertical_distance >= best_cost) {
                        break;
                    }
                    
                    // Try sub-rows in this row
                    auto [sub_row, result] = find_best_sub_row_in_row(
                        cell_idx, row_idx, current_data, best_cost, add_penalty);
                    
                    if (result.valid && result.cost < best_cost) {
                        best_cost = result.cost;
                        best_sub_row = sub_row;
                        best_result = result;
                        found_valid_placement = true;
                    }
                }
                
                // If found valid placement with penalty, no need to try without
                if (found_valid_placement) {
                    break;
                }
            }
            
            // Place cell in best sub-row (if found)
            if (found_valid_placement && best_sub_row) {
                best_sub_row->add_cell(cell_idx, current_data.cells);
                
                // Use incremental final placement
                PlaceRow::place_row_incremental_final(best_sub_row, cell_idx, 
                                                     current_data.cells, best_result);
                current_data.cells[cell_idx].is_legalized = true;
            } else {
                success = false;
                break;
            }
        }
        
        if (success) {
            std::cout << "\nApplying site alignment..." << std::endl;
            apply_site_alignment_all(current_data);
            
            double total_disp = current_data.calculate_total_displacement();
            
            if (total_disp < best_total_displacement) {
                best_total_displacement = total_disp;
                best_data = current_data;
                found_valid_solution = true;
            }
        }
    }
    
    if (!found_valid_solution) {
        return false;
    }
    
    data = best_data;
    return true;
}

std::pair<SubRow*, PlacementResult> AbacusLegalizer::find_best_sub_row_in_row(
    int cell_idx, int row_idx, const PlacementData& data, 
    double current_best_cost, bool add_penalty) {
    
    SubRow* best_sub_row = nullptr;
    PlacementResult best_result;
    best_result.valid = false;
    best_result.cost = std::numeric_limits<double>::infinity();
    
    // Get sub-rows for this row
    std::vector<SubRow*> sub_rows = get_sub_rows_for_row(data, row_idx);
    
    // Sort sub-rows by their distance to cell's x position for better pruning
    const Cell& cell = data.cells[cell_idx];
    std::sort(sub_rows.begin(), sub_rows.end(), 
        [&cell](SubRow* a, SubRow* b) {
            double dist_a = 0;
            if (cell.original_x < a->start_x) {
                dist_a = a->start_x - cell.original_x;
            } else if (cell.original_x > a->end_x) {
                dist_a = cell.original_x - a->end_x;
            }
            
            double dist_b = 0;
            if (cell.original_x < b->start_x) {
                dist_b = b->start_x - cell.original_x;
            } else if (cell.original_x > b->end_x) {
                dist_b = cell.original_x - b->end_x;
            }
            
            return dist_a < dist_b;
        });
    
    // Try each sub-row
    for (SubRow* sub_row : sub_rows) {
        PlacementResult result = try_place_cell_trial(cell_idx, sub_row, data, 
                                                     current_best_cost, add_penalty);
        
        if (result.valid && result.cost < best_result.cost) {
            best_result = result;
            best_sub_row = sub_row;
        }
    }
    
    return {best_sub_row, best_result};
}

// Keep other methods unchanged
std::vector<SubRow*> AbacusLegalizer::get_sub_rows_for_row(const PlacementData& data, int row_idx) {
    std::vector<SubRow*> sub_rows;
    if (row_idx >= 0 && row_idx < static_cast<int>(data.rows.size())) {
        for (auto& sub_row : data.rows[row_idx].sub_rows) {
            sub_rows.push_back(const_cast<SubRow*>(&sub_row));
        }
    }
    return sub_rows;
}

int AbacusLegalizer::find_nearest_row(const Cell& cell, const PlacementData& data) {
    int nearest_row = 0;
    double min_distance = std::numeric_limits<double>::infinity();
    
    for (int i = 0; i < static_cast<int>(data.rows.size()); ++i) {
        double distance = std::abs(data.rows[i].y - cell.original_y);
        if (distance < min_distance) {
            min_distance = distance;
            nearest_row = i;
        }
    }
    
    return nearest_row;
}

void AbacusLegalizer::apply_site_alignment_all(PlacementData& data) {
    // Apply site alignment to all cells after placement
    for (auto& row : data.rows) {
        for (auto& sub_row : row.sub_rows) {
            if (sub_row.cells.empty()) continue;
            
            // Sort cells by current x position
            std::vector<int> cell_indices = sub_row.cells;
            std::sort(cell_indices.begin(), cell_indices.end(),
                      [&data](int a, int b) {
                          return data.cells[a].current_x < data.cells[b].current_x;
                      });
            
            // Apply site alignment
            for (int i = 0; i < static_cast<int>(cell_indices.size()); ++i) {
                int cell_idx = cell_indices[i];
                Cell& cell = data.cells[cell_idx];
                
                // Align to site boundary  
                double relative_x = cell.current_x - sub_row.start_x;
                int site_offset = static_cast<int>(std::round(relative_x / sub_row.site_width));
                double new_x = sub_row.start_x + site_offset * sub_row.site_width;
                
                cell.current_x = new_x;
            }
            
            // Fix any overlaps caused by site alignment
            for (int i = 1; i < static_cast<int>(cell_indices.size()); ++i) {
                Cell& prev_cell = data.cells[cell_indices[i-1]];
                Cell& curr_cell = data.cells[cell_indices[i]];
                
                double prev_end = prev_cell.current_x + prev_cell.width;
                if (prev_end > curr_cell.current_x + PlaceRow::EPSILON) {
                    // Move current cell to next available site
                    double required_start = prev_end;
                    double relative_x = required_start - sub_row.start_x;
                    int site_offset = static_cast<int>(std::ceil(relative_x / sub_row.site_width));
                    double new_x = sub_row.start_x + site_offset * sub_row.site_width;
                    
                    curr_cell.current_x = new_x;
                }
            }
        }
    }
}

int AbacusLegalizer::get_sub_row_index(SubRow* sub_row, const PlacementData& data) {
    for (int i = 0; i < static_cast<int>(data.all_sub_rows.size()); ++i) {
        if (data.all_sub_rows[i] == sub_row) {
            return i;
        }
    }
    return -1;
}

double AbacusLegalizer::calculate_lower_bound_cost(const Cell& cell, const SubRow& sub_row) {
    return std::abs(sub_row.y - cell.original_y);
}

std::vector<int> AbacusLegalizer::sort_cells_by_x(const std::vector<Cell>& cells, bool ascending) {
    std::vector<int> indices(cells.size());
    for (int i = 0; i < static_cast<int>(cells.size()); ++i) {
        indices[i] = i;
    }
    
    if (!ascending) {
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