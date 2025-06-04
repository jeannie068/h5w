#include "place_row.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>

PlacementResult PlaceRow::place_row_trial(SubRow* sub_row, int new_cell_index,
                                         const std::vector<Cell>& cells, 
                                         double max_displacement_constraint) {
    PlacementResult result;
    result.valid = false;
    result.cost = std::numeric_limits<double>::infinity();
    
    if (!sub_row) {
        // std::cout << "      PlaceRow trial: sub_row is null" << std::endl;
        return result;
    }
    
    // std::cout << "      PlaceRow trial: cell " << cells[new_cell_index].name 
    //          << " to sub-row with " << sub_row->cells.size() << " existing cells" << std::endl;
    
    // Get all cells in this sub-row including the new one
    std::vector<int> cell_indices = sub_row->cells;
    
    // Add new cell to the list for trial
    cell_indices.push_back(new_cell_index);
    
    // Sort cells by their original x position
    std::sort(cell_indices.begin(), cell_indices.end(),
              [&cells](int a, int b) {
                  return cells[a].original_x < cells[b].original_x;
              });
    
    // std::cout << "      Sorted cell order: ";
    for (int idx : cell_indices) {
    //     std::cout << cells[idx].name << "(" << cells[idx].original_x << ") ";
    }
    // std::cout << std::endl;
    
    // Create clusters using dynamic programming approach
    std::vector<Cluster> clusters;
    
    for (int i = 0; i < static_cast<int>(cell_indices.size()); ++i) {
        int cell_index = cell_indices[i];
        
        // Adjust cell position to be within sub-row bounds (like reference implementation)
        double cell_x = cells[cell_index].original_x;
        if (cell_x < sub_row->start_x) {
            cell_x = sub_row->start_x;
        } else if (cell_x > sub_row->end_x - cells[cell_index].width) {
            cell_x = sub_row->end_x - cells[cell_index].width;
        }
        // std::cout << "      Cell " << cells[cell_index].name 
        //          << " adjusted position: " << cells[cell_index].original_x 
        //          << " -> " << cell_x << std::endl;
        
        bool create_new_cluster = clusters.empty();
        
        if (!create_new_cluster) {
            Cluster& last_cluster = clusters.back();
            double last_cluster_end = last_cluster.optimal_x + last_cluster.total_width;
            
            if (last_cluster_end <= cell_x + EPSILON) {
                create_new_cluster = true;
            }
        }
        
        if (create_new_cluster) {
            Cluster new_cluster;
            new_cluster.first_cell_index = i;
            new_cluster.last_cell_index = i;
            new_cluster.total_weight = cells[cell_index].weight;
            new_cluster.total_width = cells[cell_index].width;
            new_cluster.q_value = cells[cell_index].weight * cell_x;
            new_cluster.optimal_x = cell_x;
            clusters.push_back(new_cluster);
            // std::cout << "      Created new cluster for cell " << cells[cell_index].name 
            //          << " at x=" << cell_x << std::endl;
        } else {
            // Add cell to last cluster
            Cluster& cluster = clusters.back();
            cluster.last_cell_index = i;
            cluster.total_weight += cells[cell_index].weight;
            cluster.q_value += cells[cell_index].weight * (cell_x - cluster.total_width);
            cluster.total_width += cells[cell_index].width;
            // std::cout << "      Added cell " << cells[cell_index].name << " to existing cluster" << std::endl;
            
            // Collapse clusters if necessary
            collapse_clusters_trial(clusters, clusters.size() - 1, sub_row->start_x, sub_row->end_x);
        }
    }
    
    // std::cout << "      Total clusters: " << clusters.size() << std::endl;
    
    // Calculate positions from clusters (without site alignment for now)
    result.cell_positions.clear();
    for (const auto& cluster : clusters) {
        double current_x = cluster.optimal_x;
        // std::cout << "      Cluster starting at x=" << current_x << ", width=" << cluster.total_width << std::endl;
        
        for (int i = cluster.first_cell_index; i <= cluster.last_cell_index; ++i) {
            int cell_idx = cell_indices[i];
            result.cell_positions.push_back({cell_idx, current_x});
            // std::cout << "        Cell " << cells[cell_idx].name << " at x=" << current_x << std::endl;
            current_x += cells[cell_idx].width;
        }
    }
    
    // Check constraints
    result.valid = true;
    
    // Check boundary constraints
    for (const auto& pos : result.cell_positions) {
        double cell_end = pos.second + cells[pos.first].width;
        if (pos.second < sub_row->start_x - EPSILON || cell_end > sub_row->end_x + EPSILON) {
            // std::cout << "      Boundary constraint violated for cell " << cells[pos.first].name 
            //          << ": x=" << pos.second << ", end=" << cell_end 
            //          << ", sub_row=[" << sub_row->start_x << ", " << sub_row->end_x << "]" << std::endl;
            result.valid = false;
            return result;
        }
    }
    
    // Check max displacement constraint
    for (const auto& pos : result.cell_positions) {
        int cell_idx = pos.first;
        double dx = pos.second - cells[cell_idx].original_x;
        double dy = sub_row->y - cells[cell_idx].original_y;
        double displacement = std::sqrt(dx * dx + dy * dy);
        
        if (displacement > max_displacement_constraint) {
            // std::cout << "      Max displacement violated for cell " << cells[cell_idx].name 
            //          << ": " << displacement << " > " << max_displacement_constraint << std::endl;
            result.valid = false;
            return result;
        }
        
        // Calculate cost for the new cell
        if (cell_idx == new_cell_index) {
            result.cost = displacement;
            // std::cout << "      New cell displacement cost: " << result.cost << std::endl;
        }
    }
    
    // std::cout << "      PlaceRow trial result: valid=" << result.valid << ", cost=" << result.cost << std::endl;
    return result;
}

void PlaceRow::place_row_final(SubRow* sub_row, const PlacementResult& result,
                              std::vector<Cell>& cells) {
    if (!result.valid || !sub_row) {
        return;
    }
    
    // Apply the calculated positions
    for (const auto& pos : result.cell_positions) {
        cells[pos.first].current_x = pos.second;
        cells[pos.first].current_y = sub_row->y;
    }
}

// Helper function for trial mode cluster collapse
void PlaceRow::collapse_clusters_trial(std::vector<Cluster>& clusters, int cluster_index,
                                      double sub_row_start_x, double sub_row_end_x) {
    if (cluster_index < 0 || cluster_index >= static_cast<int>(clusters.size())) {
        return;
    }
    
    Cluster& current_cluster = clusters[cluster_index];
    
    // Calculate optimal position
    current_cluster.optimal_x = current_cluster.q_value / current_cluster.total_weight;
    // std::cout << "        Cluster optimal_x before bounds: " << current_cluster.optimal_x << std::endl;
    
    // Limit position between sub-row boundaries
    current_cluster.optimal_x = std::max(current_cluster.optimal_x, sub_row_start_x);
    current_cluster.optimal_x = std::min(current_cluster.optimal_x, 
                                       sub_row_end_x - current_cluster.total_width);
    // std::cout << "        Cluster optimal_x after bounds: " << current_cluster.optimal_x 
    //          << ", width=" << current_cluster.total_width << std::endl;
    
    // Check overlap with predecessor cluster
    if (cluster_index > 0) {
        Cluster& prev_cluster = clusters[cluster_index - 1];
        double prev_cluster_end = prev_cluster.optimal_x + prev_cluster.total_width;
        
        if (prev_cluster_end > current_cluster.optimal_x + EPSILON) {
            // std::cout << "        Clusters overlap, merging..." << std::endl;
            // Merge clusters
            prev_cluster.last_cell_index = current_cluster.last_cell_index;
            prev_cluster.total_weight += current_cluster.total_weight;
            prev_cluster.q_value += current_cluster.q_value - current_cluster.total_weight * prev_cluster.total_width;
            prev_cluster.total_width += current_cluster.total_width;
            
            clusters.erase(clusters.begin() + cluster_index);
            
            // Recursively collapse the merged cluster
            collapse_clusters_trial(clusters, cluster_index - 1, sub_row_start_x, sub_row_end_x);
        }
    }
}

// Legacy implementation - DEPRECATED - Use place_row_trial/place_row_final instead
bool PlaceRow::place_row(SubRow* sub_row, std::vector<Cell>& cells, 
                        double max_displacement_constraint, bool trial_mode) {
    if (!sub_row) {
        return true;
    }
    
    // If no cells to place, return success
    if (sub_row->cells.empty()) {
        return true;
    }
    
    // For backward compatibility, use the new trial/final approach
    if (trial_mode) {
        // In trial mode, we shouldn't modify cells
        // This is a compatibility shim - the caller should use place_row_trial directly
        std::cerr << "Warning: Legacy place_row() called in trial mode. Use place_row_trial() instead." << std::endl;
        return false;
    }
    
    // In final mode, we need to place all cells in the sub-row
    // Get all cells and sort by original x
    std::vector<int> cell_indices = sub_row->cells;
    sort_cells_by_original_x(cell_indices, cells);
    
    // Use the cluster-based placement algorithm
    std::vector<Cluster> clusters;
    
    for (int i = 0; i < static_cast<int>(cell_indices.size()); ++i) {
        int cell_index = cell_indices[i];
        
        // Adjust cell position to be within sub-row bounds
        double cell_x = cells[cell_index].original_x;
        if (cell_x < sub_row->start_x) {
            cell_x = sub_row->start_x;
        } else if (cell_x > sub_row->end_x - cells[cell_index].width) {
            cell_x = sub_row->end_x - cells[cell_index].width;
        }
        
        bool create_new_cluster = clusters.empty();
        
        if (!create_new_cluster) {
            Cluster& last_cluster = clusters.back();
            double last_cluster_end = last_cluster.optimal_x + last_cluster.total_width;
            
            if (last_cluster_end <= cell_x + EPSILON) {
                create_new_cluster = true;
            }
        }
        
        if (create_new_cluster) {
            Cluster new_cluster;
            new_cluster.first_cell_index = i;
            new_cluster.last_cell_index = i;
            new_cluster.total_weight = cells[cell_index].weight;
            new_cluster.total_width = cells[cell_index].width;
            new_cluster.q_value = cells[cell_index].weight * cell_x;
            new_cluster.optimal_x = cell_x;
            clusters.push_back(new_cluster);
        } else {
            Cluster& cluster = clusters.back();
            cluster.last_cell_index = i;
            cluster.total_weight += cells[cell_index].weight;
            cluster.q_value += cells[cell_index].weight * (cell_x - cluster.total_width);
            cluster.total_width += cells[cell_index].width;
            
            collapse_clusters_trial(clusters, clusters.size() - 1, sub_row->start_x, sub_row->end_x);
        }
    }
    
    // Transform clusters to cell positions
    for (const auto& cluster : clusters) {
        double current_x = cluster.optimal_x;
        for (int i = cluster.first_cell_index; i <= cluster.last_cell_index; ++i) {
            int cell_idx = cell_indices[i];
            cells[cell_idx].current_x = current_x;
            cells[cell_idx].current_y = sub_row->y;
            current_x += cells[cell_idx].width;
        }
    }
    
    // Check constraints
    for (int cell_idx : cell_indices) {
        // Check max displacement
        double dx = cells[cell_idx].current_x - cells[cell_idx].original_x;
        double dy = cells[cell_idx].current_y - cells[cell_idx].original_y;
        double displacement = std::sqrt(dx * dx + dy * dy);
        
        if (displacement > max_displacement_constraint) {
            return false;
        }
        
        // Check boundaries
        if (cells[cell_idx].current_x < sub_row->start_x - EPSILON ||
            cells[cell_idx].current_x + cells[cell_idx].width > sub_row->end_x + EPSILON) {
            return false;
        }
    }
    
    return true;
}

// Helper implementations from original code
void PlaceRow::add_cell_to_cluster(Cluster& cluster, int cell_index, 
                                  const std::vector<Cell>& cells) {
    cluster.total_weight += cells[cell_index].weight;
    cluster.q_value += cells[cell_index].weight * (cells[cell_index].original_x - cluster.total_width);
    cluster.total_width += cells[cell_index].width;
}

void PlaceRow::add_cluster_to_cluster(Cluster& target_cluster, const Cluster& source_cluster) {
    target_cluster.last_cell_index = source_cluster.last_cell_index;
    target_cluster.total_weight += source_cluster.total_weight;
    target_cluster.q_value += source_cluster.q_value - source_cluster.total_weight * target_cluster.total_width;
    target_cluster.total_width += source_cluster.total_width;
}

void PlaceRow::sort_cells_by_original_x(std::vector<int>& cell_indices, 
                                       const std::vector<Cell>& cells) {
    std::sort(cell_indices.begin(), cell_indices.end(),
              [&cells](int a, int b) {
                  return cells[a].original_x < cells[b].original_x;
              });
}