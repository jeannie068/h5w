// abacus_legalizer.cpp - Modified to use incremental PlaceRow
#include "abacus_legalizer.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>

// Modified abacus_legalizer.cpp - Key optimized functions
bool AbacusLegalizer::legalize(PlacementData& data) {
    std::cout << "\n===== Starting Abacus Legalization =====" << std::endl;
    
    // Try both ascending and descending order
    PlacementData best_data = data;
    bool found_valid_solution = false;
    double best_total_displacement = std::numeric_limits<double>::infinity();

    // For large datasets, only try ascending order
    int max_directions = (data.cells.size() >= 220000) ? 1 : 2;
    
    for (int sort_dir = 0; sort_dir < max_directions; ++sort_dir) {
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
            Cell* cell = &current_data.cells[cell_idx];
            
            int base_row_idx = find_nearest_row(*cell, current_data);
            int best_row_idx = -1;
            int best_sub_row_idx = -1;
            
            // Try with penalty first, then without if needed
            for (int add_penalty = 1; add_penalty >= 0; --add_penalty) {
                double best_cost = std::numeric_limits<double>::max();
                
                // Search upward from base row
                for (int row_idx = base_row_idx; row_idx >= 0; --row_idx) {
                    double y_distance = std::abs(cell->original_y - current_data.rows[row_idx].y);
                    if (y_distance >= best_cost) {
                        break;
                    }
                    
                    auto [sub_row_idx, cost] = try_place_in_row(cell, row_idx, current_data, add_penalty);
                    
                    if (sub_row_idx >= 0 && cost < best_cost) {
                        best_cost = cost;
                        best_row_idx = row_idx;
                        best_sub_row_idx = sub_row_idx;
                    }
                }
                
                // Search downward from base row
                for (int row_idx = base_row_idx + 1; row_idx < static_cast<int>(current_data.rows.size()); ++row_idx) {
                    double y_distance = std::abs(cell->original_y - current_data.rows[row_idx].y);
                    if (y_distance >= best_cost) {
                        break;
                    }
                    
                    auto [sub_row_idx, cost] = try_place_in_row(cell, row_idx, current_data, add_penalty);
                    
                    if (sub_row_idx >= 0 && cost < best_cost) {
                        best_cost = cost;
                        best_row_idx = row_idx;
                        best_sub_row_idx = sub_row_idx;
                    }
                }
                
                // If found valid placement, break
                if (best_sub_row_idx >= 0) {
                    break;
                }
            }
            
            // Place cell in best sub-row
            if (best_sub_row_idx >= 0) {
                SubRow* best_sub_row = &current_data.rows[best_row_idx].sub_rows[best_sub_row_idx];
                best_sub_row->add_cell(cell_idx, current_data.cells);
                PlaceRow::place_row_final(best_sub_row, cell, current_data.cells);
                cell->is_legalized = true;
                processed++;
                
                // Progress report
                if (processed % 10000 == 0 && data.cells.size() > 50000) {
                    std::cout << "  Processed " << processed << "/" << data.cells.size() 
                              << " cells..." << std::endl;
                }
            } else {
                std::cout << "  Failed to place cell " << cell->name << std::endl;
                success = false;
                break;
            }
        }
        
        if (success) {
            // Apply final site alignment
            determine_final_positions(current_data);
            
            double total_disp = current_data.calculate_total_displacement();
            std::cout << "  Total displacement: " << total_disp << std::endl;
            
            if (total_disp < best_total_displacement) {
                best_total_displacement = total_disp;
                best_data = current_data;
                found_valid_solution = true;
            }
        }
    }
    
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

std::pair<int, double> AbacusLegalizer::try_place_in_row(Cell* cell, int row_idx, 
                                                         const PlacementData& data, 
                                                         bool add_penalty) {
    if (row_idx < 0 || row_idx >= static_cast<int>(data.rows.size())) {
        return {-1, std::numeric_limits<double>::max()};
    }
    
    const Row& row = data.rows[row_idx];
    
    // Find sub-row with minimum displacement
    int best_sub_row_idx = -1;
    double best_cost = std::numeric_limits<double>::max();
    
    for (int i = 0; i < static_cast<int>(row.sub_rows.size()); ++i) {
        const SubRow& sub_row = row.sub_rows[i];
        
        // Quick check: can cell fit?
        if (cell->width > sub_row.free_width) {
            continue;
        }
        
        // Calculate x displacement for this sub-row
        double x_displacement = 0;
        if (cell->original_x < sub_row.start_x) {
            x_displacement = sub_row.start_x - cell->original_x;
        } else if (cell->original_x + cell->width > sub_row.end_x) {
            x_displacement = cell->original_x + cell->width - sub_row.end_x;
        }
        
        // Early skip if x displacement alone is too large
        if (x_displacement > best_cost) {
            continue;
        }
        
        // Try placing in this sub-row
        auto [valid, cost] = PlaceRow::place_row_trial(
            const_cast<SubRow*>(&sub_row), cell, data.cells, 
            data.max_displacement_constraint, add_penalty);
        
        if (valid >= 0 && cost < best_cost) {
            best_cost = cost;
            best_sub_row_idx = i;
        }
    }
    
    return {best_sub_row_idx, best_cost};
}