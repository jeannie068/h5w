#include "place_row.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>

bool PlaceRow::place_row(SubRow* sub_row, std::vector<Cell>& cells, 
                        double max_displacement_constraint, bool trial_mode) {
    if (!sub_row || sub_row->cells.empty()) {
        return true;
    }
    
    // Step 1: Sort cells by their original x position (Algorithm 2 assumption)
    std::vector<int> cell_indices = sub_row->cells;
    sort_cells_by_original_x(cell_indices, cells);
    
    // Step 2: Create clusters using dynamic programming approach
    std::vector<Cluster> clusters;
    
    for (int i = 0; i < static_cast<int>(cell_indices.size()); ++i) {
        int cell_index = cell_indices[i];
        
        // Check if this is first cell or if it doesn't overlap with last cluster
        bool create_new_cluster = clusters.empty();
        
        if (!create_new_cluster) {
            Cluster& last_cluster = clusters.back();
            double last_cluster_end = calculate_cluster_optimal_x(last_cluster, cells) + 
                                    calculate_total_cluster_width(last_cluster);
            
            // If cell doesn't overlap with last cluster, create new cluster
            // Use epsilon for floating point comparison
            if (last_cluster_end <= cells[cell_index].original_x + EPSILON) {
                create_new_cluster = true;
            }
        }
        
        if (create_new_cluster) {
            // Create new cluster (Algorithm 2, lines 4-8)
            Cluster new_cluster;
            new_cluster.first_cell_index = i;  // Index in cell_indices array
            new_cluster.last_cell_index = i;   // Index in cell_indices array
            new_cluster.total_weight = 0.0;
            new_cluster.total_width = 0.0;
            new_cluster.q_value = 0.0;
            new_cluster.optimal_x = cells[cell_index].original_x;
            
            add_cell_to_cluster(new_cluster, cell_index, cells);
            clusters.push_back(new_cluster);
        } else {
            // Add cell to last cluster (Algorithm 2, lines 10-11)
            clusters.back().last_cell_index = i;  // Update last cell index in array
            add_cell_to_cluster(clusters.back(), cell_index, cells);
            
            // Collapse clusters if necessary (Algorithm 2, line 11)
            int last_cluster_index = static_cast<int>(clusters.size()) - 1;
            if (!collapse_cluster(clusters, last_cluster_index, 
                                sub_row->start_x, sub_row->end_x, sub_row->site_width)) {
                return false; // Failed to resolve overlaps
            }
        }
    }
    
    // Step 3: Apply site alignment to clusters
    if (!apply_site_alignment(clusters, sub_row)) {
        return false;
    }
    
    // Step 4: Check constraints using computed positions
    // Create temporary positions for constraint checking
    std::vector<std::pair<double, double>> temp_positions;
    std::vector<int> cell_indices_sorted = sub_row->cells;
    sort_cells_by_original_x(cell_indices_sorted, cells);
    
    // Calculate positions from clusters for constraint checking
    temp_positions.resize(cell_indices_sorted.size());
    for (const auto& cluster : clusters) {
        double current_x = cluster.optimal_x;
        for (int i = cluster.first_cell_index; i <= cluster.last_cell_index; ++i) {
            temp_positions[i] = {current_x, static_cast<double>(sub_row->y)};
            current_x += cells[cell_indices_sorted[i]].width;
        }
    }
    
    if (!check_row_boundary_constraint(clusters, sub_row->start_x, sub_row->end_x)) {
        return false;
    }
    
    if (!check_max_displacement_constraint_with_positions(cell_indices_sorted, cells, 
                                                        temp_positions, max_displacement_constraint)) {
        return false;
    }
    
    // Step 5: Transform cluster positions to cell positions (Algorithm 2, lines 14-21)
    if (!trial_mode) {
        transform_clusters_to_cell_positions(clusters, sub_row, cells);
    }
    
    return true;
}

void PlaceRow::add_cell_to_cluster(Cluster& cluster, int cell_index, 
                                  const std::vector<Cell>& cells) {
    // Update cluster properties (similar to Algorithm 2, lines 24-26)
    cluster.total_weight += cells[cell_index].weight;
    cluster.q_value += cells[cell_index].weight * (cells[cell_index].original_x - cluster.total_width);
    cluster.total_width += cells[cell_index].width;
}

void PlaceRow::add_cluster_to_cluster(Cluster& target_cluster, const Cluster& source_cluster) {
    // Merge two clusters (Algorithm 2, lines 28-31)
    target_cluster.last_cell_index = source_cluster.last_cell_index;
    target_cluster.total_weight += source_cluster.total_weight;
    target_cluster.q_value += source_cluster.q_value - source_cluster.total_weight * target_cluster.total_width;
    target_cluster.total_width += source_cluster.total_width;
}

bool PlaceRow::collapse_cluster(std::vector<Cluster>& clusters, int cluster_index,
                               double sub_row_start_x, double sub_row_end_x, int site_width) {
    if (cluster_index < 0 || cluster_index >= static_cast<int>(clusters.size())) {
        return true;
    }
    
    Cluster& current_cluster = clusters[cluster_index];
    
    // Calculate optimal position (Algorithm 2, line 33)
    current_cluster.optimal_x = current_cluster.q_value / current_cluster.total_weight;
    
    // Limit position between sub-row boundaries (Algorithm 2, lines 34-35)
    current_cluster.optimal_x = std::max(current_cluster.optimal_x, sub_row_start_x);
    current_cluster.optimal_x = std::min(current_cluster.optimal_x, 
                                       sub_row_end_x - current_cluster.total_width);
    
    // Check overlap with predecessor cluster (Algorithm 2, lines 37-41)
    if (cluster_index > 0) {
        Cluster& prev_cluster = clusters[cluster_index - 1];
        double prev_cluster_end = prev_cluster.optimal_x + prev_cluster.total_width;
        
        // Use epsilon for floating point comparison
        if (prev_cluster_end > current_cluster.optimal_x + EPSILON) {
            // Merge current cluster to predecessor
            add_cluster_to_cluster(prev_cluster, current_cluster);
            clusters.erase(clusters.begin() + cluster_index);
            
            // Recursively collapse the merged cluster
            return collapse_cluster(clusters, cluster_index - 1, 
                                  sub_row_start_x, sub_row_end_x, site_width);
        }
    }
    
    return true;
}

double PlaceRow::calculate_cluster_optimal_x(const Cluster& cluster, 
                                           const std::vector<Cell>& cells) {
    if (cluster.total_weight == 0.0) {
        return cluster.optimal_x;
    }
    return cluster.q_value / cluster.total_weight;
}

double PlaceRow::get_site_aligned_position(double x, double sub_row_start_x, int site_width) {
    // Apply site alignment using floor (as per alg.md Phase 4)
    double relative_x = x - sub_row_start_x;
    int site_offset = static_cast<int>(std::floor(relative_x / site_width));
    return sub_row_start_x + site_offset * site_width;
}

bool PlaceRow::apply_site_alignment(std::vector<Cluster>& clusters, SubRow* sub_row) {
    // Apply site alignment to each cluster's optimal position
    for (auto& cluster : clusters) {
        double aligned_x = get_site_aligned_position(cluster.optimal_x, 
                                                   sub_row->start_x, sub_row->site_width);
        
        // Ensure aligned position doesn't move cluster out of sub-row
        if (aligned_x < sub_row->start_x - EPSILON) {
            aligned_x = sub_row->start_x;
        }
        if (aligned_x + cluster.total_width > sub_row->end_x + EPSILON) {
            // Try to move to previous site
            aligned_x = get_site_aligned_position(sub_row->end_x - cluster.total_width,
                                                sub_row->start_x, sub_row->site_width);
            if (aligned_x < sub_row->start_x - EPSILON) {
                return false; // Cannot fit in sub-row with site alignment
            }
        }
        
        cluster.optimal_x = aligned_x;
    }
    
    // Check for overlaps after site alignment and resolve if necessary
    for (int i = 1; i < static_cast<int>(clusters.size()); ++i) {
        double prev_end = clusters[i-1].optimal_x + clusters[i-1].total_width;
        // Use epsilon for floating point comparison
        if (prev_end > clusters[i].optimal_x + EPSILON) {
            // Push current cluster to next site boundary
            double min_x = prev_end;
            clusters[i].optimal_x = get_site_aligned_position(min_x, sub_row->start_x, sub_row->site_width);
            
            // If aligned position is still overlapping, push to next site
            if (clusters[i].optimal_x < min_x - EPSILON) {
                clusters[i].optimal_x += sub_row->site_width;
            }
            
            // Check if pushed cluster still fits in sub-row
            if (clusters[i].optimal_x + clusters[i].total_width > sub_row->end_x + EPSILON) {
                return false;
            }
        }
    }
    
    return true;
}

bool PlaceRow::check_max_displacement_constraint(const std::vector<int>& cell_indices,
                                               const std::vector<Cell>& cells,
                                               double max_displacement_constraint) {
    for (int cell_index : cell_indices) {
        if (cells[cell_index].get_displacement() > max_displacement_constraint) {
            return false;
        }
    }
    return true;
}

bool PlaceRow::check_max_displacement_constraint_with_positions(const std::vector<int>& cell_indices,
                                                              const std::vector<Cell>& cells,
                                                              const std::vector<std::pair<double, double>>& positions,
                                                              double max_displacement_constraint) {
    for (int i = 0; i < static_cast<int>(cell_indices.size()); ++i) {
        int cell_index = cell_indices[i];
        double new_x = positions[i].first;
        double new_y = positions[i].second;
        
        // Calculate displacement with new position
        double dx = new_x - cells[cell_index].original_x;
        double dy = new_y - cells[cell_index].original_y;
        double displacement = std::sqrt(dx * dx + dy * dy);
        
        if (displacement > max_displacement_constraint) {
            return false;
        }
    }
    return true;
}

bool PlaceRow::check_row_boundary_constraint(const std::vector<Cluster>& clusters,
                                           double sub_row_start_x, double sub_row_end_x) {
    for (const auto& cluster : clusters) {
        // Use epsilon for floating point comparison
        if (cluster.optimal_x < sub_row_start_x - EPSILON || 
            cluster.optimal_x + cluster.total_width > sub_row_end_x + EPSILON) {
            return false;
        }
    }
    return true;
}

void PlaceRow::transform_clusters_to_cell_positions(const std::vector<Cluster>& clusters,
                                                  SubRow* sub_row, std::vector<Cell>& cells) {
    // Transform cluster positions to individual cell positions (Algorithm 2, lines 14-21)
    std::vector<int> cell_indices = sub_row->cells;
    sort_cells_by_original_x(cell_indices, cells);
    
    for (const auto& cluster : clusters) {
        double current_x = cluster.optimal_x;
        
        // Iterate through cells in this cluster using indices
        for (int i = cluster.first_cell_index; i <= cluster.last_cell_index; ++i) {
            int cell_index = cell_indices[i];
            cells[cell_index].current_x = current_x;
            cells[cell_index].current_y = sub_row->y;
            current_x += cells[cell_index].width;
        }
    }
}

double PlaceRow::calculate_total_cluster_width(const Cluster& cluster) {
    return cluster.total_width;
}

bool PlaceRow::validate_non_overlap(const std::vector<int>& cell_indices,
                                   const std::vector<Cell>& cells) {
    for (int i = 1; i < static_cast<int>(cell_indices.size()); ++i) {
        int current_cell = cell_indices[i];
        int prev_cell = cell_indices[i-1];
        
        double prev_end = cells[prev_cell].current_x + cells[prev_cell].width;
        // Use epsilon for floating point comparison
        if (prev_end > cells[current_cell].current_x + EPSILON) {
            return false;
        }
    }
    return true;
}

void PlaceRow::sort_cells_by_original_x(std::vector<int>& cell_indices, 
                                       const std::vector<Cell>& cells) {
    std::sort(cell_indices.begin(), cell_indices.end(),
              [&cells](int a, int b) {
                  return cells[a].original_x < cells[b].original_x;
              });
}