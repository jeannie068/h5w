#include "abacus_legalizer.hpp"
#include <algorithm>
#include <cmath>
#include <iostream> // For debugging

bool AbacusLegalizer::legalize(PlacementData& data) {
    std::cout << "\n===== Starting Abacus Legalization =====" << std::endl;
    data.initialize_sub_row_pointers(); // Ensure sub-row pointers and cell IDs are set

    PlacementData best_data_so_far = data; // Store initial state
    double best_total_displacement = std::numeric_limits<double>::infinity();
    bool found_at_least_one_solution = false;

    int num_sort_directions = 1; // Default to ascending only
     if (data.cells.size() <= 80000) { // Try both for smaller, or always for now
        num_sort_directions = 2;
     }


    for (int sort_direction_mode = 0; sort_direction_mode < num_sort_directions; ++sort_direction_mode) {
        PlacementData current_run_data = data; // Reset to original data for each sort direction
        // Reset subrow states (clusters and cell lists) for the new run
        for(auto& row : current_run_data.rows){
            for(auto& sr : row.sub_rows){
                sr.cells.clear();
                sr.last_cluster = nullptr;
            }
        }


        bool sort_ascending = (sort_direction_mode == 0);
        std::cout << "\n--- Trying sort direction: " << (sort_ascending ? "ascending" : "descending") << " ---" << std::endl;

        std::vector<int> cell_processing_order = sort_cells_by_x(current_run_data.cells, sort_ascending);
        bool current_run_successful = true;

        for (int cell_idx_in_order : cell_processing_order) {
            const Cell& cell_to_place = current_run_data.cells[cell_idx_in_order];
            if(cell_to_place.is_legalized) continue; // Should not happen if data is reset properly

            double best_cost_for_this_cell = std::numeric_limits<double>::infinity();
            SubRow* best_sub_row_for_this_cell = nullptr;
            // PlacementTrialResult best_trial_result_for_this_cell; // Not strictly needed to store full result here

            int nearest_row_original_idx = find_nearest_row_idx(cell_to_place, current_run_data);
            bool found_placement_for_this_cell = false;

            for (int penalty_check_mode = 1; penalty_check_mode >= 0; --penalty_check_mode) {
                bool check_penalty_flag = (penalty_check_mode == 1);

                // Search upwards from nearest row (including nearest)
                for (int r_idx = nearest_row_original_idx; r_idx >= 0; --r_idx) {
                    double vertical_dist = std::abs(static_cast<double>(current_run_data.rows[r_idx].y) - cell_to_place.original_y);
                    if (vertical_dist >= best_cost_for_this_cell && best_cost_for_this_cell != std::numeric_limits<double>::infinity()) {
                        break; 
                    }
                    for (size_t sr_loop_idx = 0; sr_loop_idx < current_run_data.rows[r_idx].sub_rows.size(); ++sr_loop_idx) {
                        SubRow* current_sub_row = &current_run_data.rows[r_idx].sub_rows[sr_loop_idx];
                         if (static_cast<double>(cell_to_place.width) > (current_sub_row->end_x - current_sub_row->start_x) + PlaceRow::EPSILON) continue;


                        PlacementTrialResult trial_res = PlaceRow::place_row_trial(current_sub_row, cell_idx_in_order, current_run_data, current_run_data.max_displacement_constraint, check_penalty_flag);
                        if (trial_res.valid && trial_res.cost < best_cost_for_this_cell) {
                            best_cost_for_this_cell = trial_res.cost;
                            best_sub_row_for_this_cell = current_sub_row;
                            found_placement_for_this_cell = true;
                        }
                    }
                }

                // Search downwards from nearest row (excluding nearest if already checked)
                for (int r_idx = nearest_row_original_idx + 1; r_idx < static_cast<int>(current_run_data.rows.size()); ++r_idx) {
                    double vertical_dist = std::abs(static_cast<double>(current_run_data.rows[r_idx].y) - cell_to_place.original_y);
                     if (vertical_dist >= best_cost_for_this_cell && best_cost_for_this_cell != std::numeric_limits<double>::infinity()) {
                        break;
                    }
                     for (size_t sr_loop_idx = 0; sr_loop_idx < current_run_data.rows[r_idx].sub_rows.size(); ++sr_loop_idx) {
                        SubRow* current_sub_row = &current_run_data.rows[r_idx].sub_rows[sr_loop_idx];
                        if (static_cast<double>(cell_to_place.width) > (current_sub_row->end_x - current_sub_row->start_x) + PlaceRow::EPSILON) continue;

                        PlacementTrialResult trial_res = PlaceRow::place_row_trial(current_sub_row, cell_idx_in_order, current_run_data, current_run_data.max_displacement_constraint, check_penalty_flag);
                        if (trial_res.valid && trial_res.cost < best_cost_for_this_cell) {
                            best_cost_for_this_cell = trial_res.cost;
                            best_sub_row_for_this_cell = current_sub_row;
                            found_placement_for_this_cell = true;
                        }
                    }
                }
                if (found_placement_for_this_cell) break; // Found with this penalty mode
            }


            if (best_sub_row_for_this_cell) {
                PlaceRow::place_row_final(best_sub_row_for_this_cell, cell_idx_in_order, current_run_data);
                current_run_data.cells[cell_idx_in_order].is_legalized = true;
            } else {
                std::cerr << "Error: Could not find valid placement for cell " << cell_to_place.name << std::endl;
                current_run_successful = false;
                break; 
            }
        }

        if (current_run_successful) {
            // apply_site_alignment_to_all_subrows(current_run_data); // Done by place_row_final now
            double current_total_disp = current_run_data.calculate_total_displacement();
            std::cout << "Sort direction " << (sort_ascending ? "ascending" : "descending") 
                      << " total displacement: " << current_total_disp << std::endl;
            if (current_total_disp < best_total_displacement) {
                best_total_displacement = current_total_disp;
                best_data_so_far = current_run_data;
                found_at_least_one_solution = true;
            }
        }
    }
    
    if (!found_at_least_one_solution) {
         std::cerr << "Error: No valid solution found after all attempts." << std::endl;
        return false;
    }

    data = best_data_so_far; // Commit the best result found
    std::cout << "\n===== Abacus Legalization Finished =====" << std::endl;
    std::cout << "Best total displacement: " << data.calculate_total_displacement() << std::endl;
    std::cout << "Max displacement: " << data.calculate_max_displacement() << std::endl;

    return true;
}


std::vector<int> AbacusLegalizer::sort_cells_by_x(const std::vector<Cell>& cells, bool ascending) {
    std::vector<int> indices(cells.size());
    for (size_t i = 0; i < cells.size(); ++i) {
        indices[i] = i; // Store original index
    }

    std::sort(indices.begin(), indices.end(),
              [&cells, ascending](int a_idx, int b_idx) {
                  if (ascending) {
                      return cells[a_idx].original_x < cells[b_idx].original_x;
                  } else {
                      return cells[a_idx].original_x > cells[b_idx].original_x;
                  }
              });
    return indices;
}

int AbacusLegalizer::find_nearest_row_idx(const Cell& cell, const PlacementData& data) {
    int nearest_idx = -1;
    double min_y_dist = std::numeric_limits<double>::infinity();

    for (size_t i = 0; i < data.rows.size(); ++i) {
        double y_dist = std::abs(static_cast<double>(data.rows[i].y) - cell.original_y);
        if (y_dist < min_y_dist) {
            min_y_dist = y_dist;
            nearest_idx = i;
        }
    }
    return nearest_idx;
}


// This function might be simplified or removed if AbacusLegalizer::legalize handles subrow iteration directly.
PlacementTrialResult AbacusLegalizer::try_place_cell_in_sub_row_trial(
    const SubRow* sub_row,
    int cell_idx,
    const PlacementData& placement_data,
    bool check_penalty_for_all_cluster_cells) {
    
    // Basic check: can cell even fit?
    if (static_cast<double>(placement_data.cells[cell_idx].width) > (sub_row->end_x - sub_row->start_x) + PlaceRow::EPSILON) {
        PlacementTrialResult res;
        res.valid = false;
        return res;
    }
    
    // More sophisticated check: available sites (from original user code)
    int cell_sites_needed = static_cast<int>(std::ceil(static_cast<double>(placement_data.cells[cell_idx].width) / sub_row->site_width));
    int used_sites = 0;
    for (int existing_cell_idx : sub_row->cells) { // sub_row->cells has FINALIZED cells. Trial needs different logic.
                                                // This check is problematic for trial mode if sub_row->cells isn't temporary.
                                                // The trial logic should consider the *current* state of clusters.
        // This site check needs to be based on the current sub_row->last_cluster's total_width
        // or a similar dynamic calculation of used space.
    }
    // For now, rely on PlaceRow::place_row_trial to handle placement feasibility.
    // The site check based on sub_row->cells is not correct for trial.

    return PlaceRow::place_row_trial(sub_row, cell_idx, placement_data, 
                                     placement_data.max_displacement_constraint, 
                                     check_penalty_for_all_cluster_cells);
}


// This function is largely superseded by direct iteration in AbacusLegalizer::legalize
std::pair<SubRow*, PlacementTrialResult> AbacusLegalizer::find_best_sub_row_in_row(
    int cell_idx, 
    int row_idx, 
    PlacementData& placement_data,
    const PlacementData& const_placement_data,
    double current_best_overall_cost,
    bool check_penalty_for_all_cluster_cells) {
    
    SubRow* best_sub_row_in_this_row = nullptr;
    PlacementTrialResult best_trial_for_this_row; // Default invalid

    const Row& current_row = const_placement_data.rows[row_idx];
    const Cell& cell_to_place = const_placement_data.cells[cell_idx];

    // Sort subrows in this row by proximity to cell's original_x
    std::vector<SubRow*> sub_rows_in_row_ptrs;
    for(size_t i=0; i < placement_data.rows[row_idx].sub_rows.size(); ++i) { // Use non-const for ptr
        sub_rows_in_row_ptrs.push_back(&placement_data.rows[row_idx].sub_rows[i]);
    }

    std::sort(sub_rows_in_row_ptrs.begin(), sub_rows_in_row_ptrs.end(), 
        [&cell_to_place](const SubRow* a, const SubRow* b) {
            double dist_a = 0;
            if (cell_to_place.original_x < a->start_x) dist_a = a->start_x - cell_to_place.original_x;
            else if (cell_to_place.original_x > a->end_x - cell_to_place.width) dist_a = cell_to_place.original_x - (a->end_x - cell_to_place.width);
            
            double dist_b = 0;
            if (cell_to_place.original_x < b->start_x) dist_b = b->start_x - cell_to_place.original_x;
            else if (cell_to_place.original_x > b->end_x - cell_to_place.width) dist_b = cell_to_place.original_x - (b->end_x - cell_to_place.width);
            return dist_a < dist_b;
        });


    for (SubRow* sub_row_ptr : sub_rows_in_row_ptrs) {
        // Check if cell can physically fit (width)
        if (static_cast<double>(cell_to_place.width) > (sub_row_ptr->end_x - sub_row_ptr->start_x) + PlaceRow::EPSILON) {
            continue;
        }

        PlacementTrialResult trial_res = PlaceRow::place_row_trial(
            sub_row_ptr, cell_idx, const_placement_data, 
            const_placement_data.max_displacement_constraint, 
            check_penalty_for_all_cluster_cells);

        if (trial_res.valid && trial_res.cost < best_trial_for_this_row.cost) {
            best_trial_for_this_row = trial_res;
            best_sub_row_in_this_row = sub_row_ptr; // Store pointer to the subrow in non-const data
        }
    }
    return {best_sub_row_in_this_row, best_trial_for_this_row};
}

void AbacusLegalizer::apply_site_alignment_to_all_subrows(PlacementData& data) {
    for (Row& row : data.rows) {
        for (SubRow& sub_row : row.sub_rows) {
            // This function is now effectively handled by PlaceRow::determine_and_apply_final_positions
            // called after each PlaceRow::place_row_final.
            // If a final sweep is desired, it could be implemented here, but might be redundant.
            // For now, this function can be left empty or removed if confidence in incremental alignment is high.
        }
    }
}
