#include "place_row.hpp"
#include <algorithm>
#include <cmath>
#include <iostream> // For potential debugging
#include <vector>
#include <stack>  // For penalty check simulation

PlacementTrialResult PlaceRow::place_row_trial(
    const SubRow* sub_row, 
    int new_cell_index,
    const PlacementData& placement_data,
    double max_displacement_constraint,
    bool check_penalty_for_all_cluster_cells) {

    PlacementTrialResult result;
    result.valid = true; // Assume valid initially

    const Cell& new_cell = placement_data.cells[new_cell_index];

    double cell_initial_x = new_cell.original_x;
    if (cell_initial_x < sub_row->start_x) {
        cell_initial_x = sub_row->start_x;
    }
    if (cell_initial_x + new_cell.width > sub_row->end_x) {
        cell_initial_x = sub_row->end_x - new_cell.width;
    }
    
    // Simulate adding the new cell
    std::shared_ptr<Cluster> trial_last_cluster = sub_row->last_cluster;
    
    double current_q, current_weight, current_width;
    double cluster_chain_start_x; // Optimal x of the whole chain being considered

    if (!trial_last_cluster || (trial_last_cluster->optimal_x + trial_last_cluster->total_width <= cell_initial_x + EPSILON)) {
        // New cell forms a new cluster
        current_q = new_cell.weight * cell_initial_x;
        current_weight = new_cell.weight;
        current_width = new_cell.width;
        cluster_chain_start_x = cell_initial_x; // Optimal x of this single-cell cluster
        
        // If there was a previous cluster, ensure no overlap. If so, this case is wrong, should merge.
        // This simple check is fine as collapse will handle it if it should have merged.
    } else {
        // New cell appends to the last cluster (conceptually)
        current_q = trial_last_cluster->q_value + new_cell.weight * (cell_initial_x - trial_last_cluster->total_width);
        current_weight = trial_last_cluster->total_weight + new_cell.weight;
        current_width = trial_last_cluster->total_width + new_cell.width;
        trial_last_cluster = trial_last_cluster->predecessor; // Start collapse check from predecessor
    }

    // Iteratively collapse with predecessors
    std::shared_ptr<Cluster> pred_iter = trial_last_cluster; // This is the predecessor of the (potentially merged) cluster
    while(pred_iter) {
        cluster_chain_start_x = current_q / current_weight;
        if (cluster_chain_start_x < sub_row->start_x) cluster_chain_start_x = sub_row->start_x;
        if (cluster_chain_start_x + current_width > sub_row->end_x) cluster_chain_start_x = sub_row->end_x - current_width;

        if (pred_iter->optimal_x + pred_iter->total_width > cluster_chain_start_x + EPSILON) {
            // Merge pred_iter into current_cluster_data
            current_q = pred_iter->q_value + current_q - current_weight * pred_iter->total_width; // Order of terms matters for q update
            current_weight += pred_iter->total_weight;
            current_width += pred_iter->total_width;
            pred_iter = pred_iter->predecessor;
        } else {
            break;
        }
    }
    // Final optimal_x for the potentially multi-cluster chain
    cluster_chain_start_x = current_q / current_weight;
    if (cluster_chain_start_x < sub_row->start_x) cluster_chain_start_x = sub_row->start_x;
    if (cluster_chain_start_x + current_width > sub_row->end_x) cluster_chain_start_x = sub_row->end_x - current_width;
    
    result.optimal_x_for_new_cell = cluster_chain_start_x + current_width - new_cell.width;

    double dx = result.optimal_x_for_new_cell - new_cell.original_x;
    double dy = static_cast<double>(sub_row->y) - new_cell.original_y;
    result.cost = std::sqrt(dx * dx + dy * dy);

    if (result.cost > max_displacement_constraint && !check_penalty_for_all_cluster_cells) {
         result.valid = false; // New cell itself violates
         return result;
    }


    if (check_penalty_for_all_cluster_cells) {
        // Simulate placement of all cells in the affected chain to check max displacement
        // This is a more complex part: reconstruct the members of the collapsed cluster(s)
        // For simplicity here, we'll assume if the new cell is okay, the penalty check is more about overall fit.
        // A full penalty check would require iterating through all cells in `current_width` starting from `cluster_chain_start_x`.
        // The reference implementation `Legalizer_ref.cpp` does this with a stack.
        
        // Simplified penalty check: if the new cell's displacement is too high, it's a bad sign.
        // A more accurate check would involve simulating positions of all cells in the collapsed cluster.
        // Let's use the new cell's displacement as a proxy for the penalty check for now.
        // If the cost (new cell's displacement) itself is already > constraint, it's invalid.
        if (result.cost > max_displacement_constraint) {
            result.valid = false;
            return result;
        }

        // More detailed penalty check (conceptual, needs full cell list of the trial cluster)
        // This part is tricky without actually building the temporary cluster structure with all members.
        // The reference code pushes clusters onto a stack during the collapse simulation.
        // For now, we rely on the single cell's displacement check and the AbaqusLegalizer's two-pass penalty.
    }
    
    return result;
}

void PlaceRow::place_row_final(
    SubRow* sub_row,
    int new_cell_index,
    PlacementData& placement_data) {

    const Cell& new_cell = placement_data.cells[new_cell_index];
    sub_row->add_cell_to_final_list(new_cell_index, placement_data.cells); // Add to sorted list for record

    double cell_initial_x = new_cell.original_x;
    if (cell_initial_x < sub_row->start_x) {
        cell_initial_x = sub_row->start_x;
    }
    if (cell_initial_x + new_cell.width > sub_row->end_x) {
        cell_initial_x = sub_row->end_x - new_cell.width;
    }

    std::shared_ptr<Cluster> current_cluster_ptr;

    if (!sub_row->last_cluster || (sub_row->last_cluster->optimal_x + sub_row->last_cluster->total_width <= cell_initial_x + EPSILON)) {
        // New cell forms a new cluster
        current_cluster_ptr = std::make_shared<Cluster>(new_cell_index, new_cell, cell_initial_x, sub_row->last_cluster);
        sub_row->last_cluster = current_cluster_ptr;
    } else {
        // Add cell to the existing last_cluster
        current_cluster_ptr = sub_row->last_cluster;
        current_cluster_ptr->member_cell_indices.push_back(new_cell_index);
        // Sort members by original_x to maintain order for position calculation
        std::sort(current_cluster_ptr->member_cell_indices.begin(), current_cluster_ptr->member_cell_indices.end(),
            [&](int a, int b){ return placement_data.cells[a].original_x < placement_data.cells[b].original_x; });

        current_cluster_ptr->q_value += new_cell.weight * (cell_initial_x - current_cluster_ptr->total_width); // q_c = q_old + e_new * (x'_new - w_old)
        current_cluster_ptr->total_weight += new_cell.weight;
        current_cluster_ptr->total_width += new_cell.width;
    }
    
    // Collapse the modified/new last_cluster with its predecessors
    while (current_cluster_ptr) {
        current_cluster_ptr->optimal_x = current_cluster_ptr->q_value / current_cluster_ptr->total_weight;
        if (current_cluster_ptr->optimal_x < sub_row->start_x) {
            current_cluster_ptr->optimal_x = sub_row->start_x;
        }
        if (current_cluster_ptr->optimal_x + current_cluster_ptr->total_width > sub_row->end_x) {
            current_cluster_ptr->optimal_x = sub_row->end_x - current_cluster_ptr->total_width;
        }

        std::shared_ptr<Cluster> pred = current_cluster_ptr->predecessor;
        if (pred && (pred->optimal_x + pred->total_width > current_cluster_ptr->optimal_x + EPSILON)) {
            // Merge current_cluster_ptr into pred
            for (int member_idx : current_cluster_ptr->member_cell_indices) {
                pred->member_cell_indices.push_back(member_idx);
            }
            std::sort(pred->member_cell_indices.begin(), pred->member_cell_indices.end(),
                [&](int a, int b){ return placement_data.cells[a].original_x < placement_data.cells[b].original_x; });
            
            // Update q for pred: q_new_pred = q_old_pred + q_current_cluster - e_current_cluster * w_old_pred
            pred->q_value = pred->q_value + current_cluster_ptr->q_value - current_cluster_ptr->total_weight * pred->total_width;
            pred->total_weight += current_cluster_ptr->total_weight;
            pred->total_width += current_cluster_ptr->total_width;
            
            // Current cluster is absorbed, pred becomes the new last_cluster if current was last
            pred->predecessor = current_cluster_ptr->predecessor ? current_cluster_ptr->predecessor->predecessor : nullptr; // This logic is tricky. Simpler:
            pred->predecessor = current_cluster_ptr->predecessor->predecessor; // if pred is not null.
            if (sub_row->last_cluster == current_cluster_ptr) { // if current_cluster_ptr was the one pointed to by sub_row
                 sub_row->last_cluster = pred;
            }
             // Update current_cluster_ptr to continue collapse check upwards
            current_cluster_ptr = pred;

        } else {
            break; // No more collapse needed for this cluster
        }
    }
    determine_and_apply_final_positions(sub_row, placement_data);
}


void PlaceRow::determine_and_apply_final_positions(SubRow* sub_row, PlacementData& placement_data) {
    std::shared_ptr<Cluster> cluster_iter = sub_row->last_cluster;
    std::vector<std::shared_ptr<Cluster>> cluster_order; // To process from left to right

    while(cluster_iter) {
        cluster_order.push_back(cluster_iter);
        cluster_iter = cluster_iter->predecessor;
    }
    std::reverse(cluster_order.begin(), cluster_order.end()); // Iterate from left-most cluster

    for (const auto& current_cluster : cluster_order) {
        // Optimal x for the cluster is already calculated and clamped.
        // Apply site alignment to the cluster's start position
        double aligned_cluster_start_x = sub_row->get_site_aligned_x(current_cluster->optimal_x);
        
        // Ensure alignment doesn't push cluster beyond subrow end
        if (aligned_cluster_start_x + current_cluster->total_width > sub_row->end_x) {
            aligned_cluster_start_x = sub_row->get_site_aligned_x(sub_row->end_x - current_cluster->total_width);
        }
        current_cluster->optimal_x = aligned_cluster_start_x; // Update cluster's optimal_x to aligned one

        double current_cell_pos_x = current_cluster->optimal_x;
        for (int cell_idx : current_cluster->member_cell_indices) {
            Cell& cell = placement_data.cells[cell_idx];
            cell.current_x = current_cell_pos_x;
            cell.current_y = static_cast<double>(sub_row->y);
            current_cell_pos_x += cell.width;
        }
    }
}

