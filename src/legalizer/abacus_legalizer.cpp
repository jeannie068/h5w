#include "abacus_legalizer.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>

bool AbacusLegalizer::legalize(PlacementData& data) {
    std::cout << "\n===== Starting Abacus Legalization =====" << std::endl;
    
    // Try both ascending and descending order at the algorithm level
    PlacementData best_data = data;
    bool found_valid_solution = false;
    double best_total_displacement = std::numeric_limits<double>::infinity();
    
    for (int sort_dir = 0; sort_dir < 1; ++sort_dir) {
        std::cout << "\n--- Trying sort direction: " << (sort_dir == 0 ? "ascending" : "descending") << " ---" << std::endl;
        
        PlacementData current_data = data;
        bool success = true;
        
        // Sort cells by x coordinate
        std::vector<int> cell_order = sort_cells_by_x(current_data.cells, sort_dir == 0);
        
        // Legalize each cell one by one
        for (int i = 0; i < static_cast<int>(cell_order.size()); ++i) {
            int cell_idx = cell_order[i];
            // std::cout << "\nProcessing cell " << current_data.cells[cell_idx].name 
            //          << " (index=" << cell_idx << ", width=" << current_data.cells[cell_idx].width
            //          << ", orig_pos=(" << current_data.cells[cell_idx].original_x 
            //          << ", " << current_data.cells[cell_idx].original_y << "))" << std::endl;
            
            double best_cost = std::numeric_limits<double>::infinity();
            int best_sub_row_idx = -1;
            PlacementResult best_result;
            
            // Try all sub-rows and find the best one
            bool found_valid_placement = false;
            
            for (int sub_row_idx = 0; sub_row_idx < static_cast<int>(current_data.all_sub_rows.size()); ++sub_row_idx) {
                SubRow* sub_row = current_data.all_sub_rows[sub_row_idx];
                
                // std::cout << "  Trying sub-row " << sub_row_idx 
                //          << " (y=" << sub_row->y 
                //          << ", x_range=[" << sub_row->start_x << ", " << sub_row->end_x << "]"
                //          << ", cells_count=" << sub_row->cells.size() << ")" << std::endl;
                
                // Check if cell can fit in sub-row
                PlacementResult result = try_place_cell_trial(cell_idx, sub_row, 
                                                             current_data, best_cost);
                
                if (result.valid) {
                    // std::cout << "    Valid placement found with cost=" << result.cost << std::endl;
                    if (result.cost < best_cost) {
                        best_cost = result.cost;
                        best_sub_row_idx = sub_row_idx;
                        best_result = result;
                        found_valid_placement = true;
                    }
                } else {
                    // std::cout << "    Invalid placement" << std::endl;
                }
            }
            
            // Place cell in best sub-row (if found)
            if (found_valid_placement && best_sub_row_idx >= 0) {
                SubRow* best_sub_row = current_data.all_sub_rows[best_sub_row_idx];
                // std::cout << "  Placing cell in sub-row " << best_sub_row_idx 
                //          << " with cost=" << best_cost << std::endl;
                
                best_sub_row->add_cell(cell_idx, current_data.cells);
                
                // Apply final placement
                PlaceRow::place_row_final(best_sub_row, best_result, current_data.cells);
                current_data.cells[cell_idx].is_legalized = true;
            } else {
                // std::cout << "  ERROR: No valid placement found for cell " 
                //          << current_data.cells[cell_idx].name << std::endl;
                success = false;
                break;
            }
        }
        
        if (success) {
            std::cout << "\nApplying site alignment..." << std::endl;
            // Apply site alignment after all cells are placed
            apply_site_alignment_all(current_data);
            
            double total_disp = current_data.calculate_total_displacement();
            // std::cout << "Total displacement for this direction: " << total_disp << std::endl;
            
            if (total_disp < best_total_displacement) {
                best_total_displacement = total_disp;
                best_data = current_data;
                found_valid_solution = true;
            }
        }
    }
    
    if (!found_valid_solution) {
        // std::cerr << "\nError: No valid placement found" << std::endl;
        return false;
    }
    
    // std::cout << "\nLegalization successful! Best total displacement: " << best_total_displacement << std::endl;
    data = best_data;
    return true;
}

PlacementResult AbacusLegalizer::try_place_cell_trial(int cell_index, SubRow* sub_row,
                                                     const PlacementData& data, 
                                                     double current_best_cost) {
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
    
    // std::cout << "    Sub-row site check: used_sites=" << used_sites 
    //          << ", available_sites=" << available_sites 
    //          << ", remaining_sites=" << remaining_sites
    //          << ", cell_sites_needed=" << cell_sites_needed << std::endl;
    
    // Check if cell can fit in terms of sites
    if (cell_sites_needed > remaining_sites) {
        // std::cout << "    Cell won't fit (site check failed)" << std::endl;
        return result;
    }
    
    // Use trial mode of PlaceRow
    result = PlaceRow::place_row_trial(sub_row, cell_index, data.cells, 
                                     data.max_displacement_constraint);
    
    if (!result.valid) {
        // std::cout << "    PlaceRow trial failed" << std::endl;
    } else {
        // std::cout << "    PlaceRow trial succeeded, cost=" << result.cost << std::endl;
    }
    
    return result;
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
                
                // std::cout << "  Aligning cell " << cell.name 
                //          << ": " << cell.current_x << " -> " << new_x << std::endl;
                
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
                    
                    // std::cout << "  Fixing overlap: moving cell " << curr_cell.name 
                    //          << " from " << curr_cell.current_x << " to " << new_x << std::endl;
                    
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