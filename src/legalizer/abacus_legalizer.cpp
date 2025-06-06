// abacus_legalizer.cpp - Modified to use incremental PlaceRow
#include "abacus_legalizer.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>

PlacementResult AbacusLegalizer::try_place_cell_trial(int cell_index, SubRow* sub_row,
                                                     const PlacementData& data, 
                                                     double current_best_cost,
                                                     bool add_penalty) {
    // Direct call to PlaceRow
    return PlaceRow::place_row_incremental_trial(sub_row, cell_index, data.cells, 
                                                data.max_displacement_constraint, add_penalty);
}

// Modified abacus_legalizer.cpp - Key optimized functions

// Fixed legalize function in abacus_legalizer.cpp

bool AbacusLegalizer::legalize(PlacementData& data) {
    std::cout << "\n===== Starting Abacus Legalization =====" << std::endl;
    
    // Try both ascending and descending order
    PlacementData best_data = data;
    bool found_valid_solution = false;
    double best_total_displacement = std::numeric_limits<double>::infinity();

    // For large datasets, only try ascending order
    // int max_dir_num = (data.cells.size() >= 220000) ? 1 : 2;
    int sort_dir = 0; // 0 for ascending, 1 for descending
    
    // for (int sort_dir = 0; sort_dir < 1; ++sort_dir) {
    std::cout << "\n--- Trying sort direction: " << (sort_dir == 0 ? "ascending" : "descending") << " ---" << std::endl;
    
    PlacementData current_data = data;
    
    // Initialize all sub-rows
    for (auto& row : current_data.rows) {
        for (auto& sub_row : row.sub_rows) {
            sub_row.last_cluster = nullptr;
            sub_row.cells.clear();
            sub_row.free_width = static_cast<int>(sub_row.end_x - sub_row.start_x);
        }
    }
    
    bool success = true;
    
    // Sort cells by x coordinate
    std::vector<int> cell_order = sort_cells_by_x(current_data.cells, sort_dir == 0);
    
    // Legalize each cell one by one
    int processed = 0;
    for (int cell_idx : cell_order) {
        double best_cost = std::numeric_limits<double>::infinity();
        SubRow* best_sub_row = nullptr;
        PlacementResult best_result;
        
        // Find nearest row based on Y coordinate
        int base_row_idx = find_nearest_row(current_data.cells[cell_idx], current_data);
        bool found_valid_placement = false;
        
        // Following reference: try with penalty first (add_penalty = 1)
        for (int add_penalty = 1; add_penalty >= 0; --add_penalty) {
            // Search from base row upward
            for (int row_idx = base_row_idx; row_idx >= 0; --row_idx) {
                // Early termination based on vertical distance
                double y_distance = std::abs(current_data.cells[cell_idx].original_y - 
                                            current_data.rows[row_idx].y);
                if (y_distance >= best_cost) {
                    break;
                }
                
                // Try all sub-rows in this row
                auto [sub_row, result] = find_best_sub_row_in_row(
                    cell_idx, row_idx, current_data, best_cost, add_penalty == 1);
                
                if (result.valid && result.cost < best_cost) {
                    best_cost = result.cost;
                    best_sub_row = sub_row;
                    best_result = result;
                    found_valid_placement = true;
                }
            }
            
            // Search from base row downward
            for (int row_idx = base_row_idx + 1; row_idx < static_cast<int>(current_data.rows.size()); ++row_idx) {
                // Early termination based on vertical distance
                double y_distance = std::abs(current_data.cells[cell_idx].original_y - 
                                            current_data.rows[row_idx].y);
                if (y_distance >= best_cost) {
                    break;
                }
                
                // Try all sub-rows in this row
                auto [sub_row, result] = find_best_sub_row_in_row(
                    cell_idx, row_idx, current_data, best_cost, add_penalty == 1);
                
                if (result.valid && result.cost < best_cost) {
                    best_cost = result.cost;
                    best_sub_row = sub_row;
                    best_result = result;
                    found_valid_placement = true;
                }
            }
            
            // If found valid placement, no need to try without penalty
            if (found_valid_placement) {
                break;
            }
        }
        
        // Place cell in best sub-row
        if (found_valid_placement && best_sub_row) {
            best_sub_row->add_cell(cell_idx, current_data.cells);
            PlaceRow::place_row_incremental_final(best_sub_row, cell_idx, 
                                                    current_data.cells, best_result);
            current_data.cells[cell_idx].is_legalized = true;
            processed++;
            
            // Progress report for large datasets
            if (processed % 100000 == 0 && data.cells.size() > 100000) {
                std::cout << "  Processed " << processed << "/" << data.cells.size() 
                            << " cells..." << std::endl;
            }
        } else {
            std::cout << "  Failed to place cell " << current_data.cells[cell_idx].name 
                        << " (index: " << cell_idx << ")" << std::endl;
            success = false;
            break;
        }
    }
    
    if (success) {
        // Apply final site alignment (following reference implementation)
        determine_final_positions(current_data);
        
        double total_disp = current_data.calculate_total_displacement();
        std::cout << "  Total displacement: " << total_disp << std::endl;
        
        if (total_disp < best_total_displacement) {
            best_total_displacement = total_disp;
            best_data = current_data;
            found_valid_solution = true;
        }
    }
    // }
    
    if (!found_valid_solution) {
        std::cout << "Legalization failed!" << std::endl;
        return false;
    }
    
    data = best_data;
    std::cout << "Best total displacement: " << best_total_displacement << std::endl;
    return true;
}

// Add new function for final position determination (like reference's determinePosition)
void AbacusLegalizer::determine_final_positions(PlacementData& data) {
    for (auto& row : data.rows) {
        for (auto& sub_row : row.sub_rows) {
            Cluster::ptr cluster = sub_row.last_cluster;
            while (cluster) {
                int x = PlaceRow::get_site_x(cluster->x, sub_row.start_x, sub_row.site_width);
                for (int cell_idx : cluster->member) {
                    data.cells[cell_idx].current_x = x;
                    data.cells[cell_idx].current_y = sub_row.y;
                    x += data.cells[cell_idx].width;
                }
                cluster = cluster->predecessor;
            }
        }
    }
}


// Fixed find_best_sub_row_in_row in abacus_legalizer.cpp

std::pair<SubRow*, PlacementResult> AbacusLegalizer::find_best_sub_row_in_row(
    int cell_idx, int row_idx, const PlacementData& data, 
    double current_best_cost, bool add_penalty) {
    
    SubRow* best_sub_row = nullptr;
    PlacementResult best_result;
    best_result.valid = false;
    best_result.cost = std::numeric_limits<double>::infinity();
    
    const Cell& cell = data.cells[cell_idx];
    
    // Get sub-rows for this row (following reference's getSubRowIdx logic)
    if (row_idx < 0 || row_idx >= static_cast<int>(data.rows.size())) {
        return {best_sub_row, best_result};
    }
    
    const Row& row = data.rows[row_idx];
    
    // Find candidate sub-rows that have enough free width
    std::vector<SubRow*> candidate_subrows;
    for (auto& sub_row : row.sub_rows) {
        if (cell.width <= sub_row.free_width) {
            candidate_subrows.push_back(const_cast<SubRow*>(&sub_row));
        }
    }
    
    if (candidate_subrows.empty()) {
        return {best_sub_row, best_result};
    }
    
    // Find sub-row with minimum x displacement (following reference logic)
    double min_x_displacement = std::numeric_limits<double>::max();
    SubRow* best_candidate = nullptr;
    
    for (SubRow* sub_row : candidate_subrows) {
        double x_displacement = 0;
        
        if (cell.original_x < sub_row->start_x) {
            x_displacement = sub_row->start_x - cell.original_x;
        } else if (cell.original_x + cell.width > sub_row->end_x) {
            x_displacement = cell.original_x + cell.width - sub_row->end_x;
        }
        // else x_displacement = 0 (cell is within sub-row range)
        
        if (x_displacement < min_x_displacement) {
            min_x_displacement = x_displacement;
            best_candidate = sub_row;
        }
    }
    
    // Try placing in the best candidate sub-row
    if (best_candidate) {
        PlacementResult result = PlaceRow::place_row_incremental_trial(
            best_candidate, cell_idx, data.cells, 
            data.max_displacement_constraint, add_penalty);
        
        if (result.valid) {
            best_result = result;
            best_sub_row = best_candidate;
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

// Fixed sort_cells_by_x function in abacus_legalizer.cpp
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