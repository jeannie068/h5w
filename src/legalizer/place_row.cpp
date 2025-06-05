// place_row.cpp - Modified key functions
#include "place_row.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stack>

Cluster::ptr PlaceRow::collapse_cluster_incremental(Cluster::ptr cluster,
                                                   double sub_row_start_x, 
                                                   double sub_row_end_x,
                                                   int site_width) {
    while (cluster) {
        // Calculate optimal position
        cluster->x = cluster->q_value / cluster->total_weight;
        
        // Limit position within sub-row boundaries
        if (cluster->x < sub_row_start_x) {
            cluster->x = sub_row_start_x;
        }
        if (cluster->x > sub_row_end_x - cluster->total_width) {
            cluster->x = sub_row_end_x - cluster->total_width;
        }
        
        // Check overlap with predecessor
        Cluster::ptr prev = cluster->predecessor;
        if (prev && prev->x + prev->total_width > cluster->x) {
            // Merge with predecessor
            prev->member.insert(prev->member.end(), cluster->member.begin(), cluster->member.end());
            prev->total_weight += cluster->total_weight;
            prev->q_value += cluster->q_value - cluster->total_weight * prev->total_width;
            prev->total_width += cluster->total_width;
            
            cluster = prev;
        } else {
            break;
        }
    }
    
    return cluster;
}

bool PlaceRow::check_cluster_constraints_with_site_alignment(
    Cluster::ptr cluster,
    const std::vector<Cell>& cells,
    double max_displacement_constraint,
    int sub_row_y,
    double sub_row_start_x,
    int site_width,
    bool add_penalty) {
    
    if (!add_penalty || !cluster) {
        return true;
    }
    
    // Collect all clusters that might have been affected
    std::vector<Cluster::ptr> affected_clusters;
    Cluster::ptr current = cluster;
    while (current) {
        affected_clusters.push_back(current);
        current = current->predecessor;
    }
    
    // Check constraints for all cells in affected clusters
    for (auto it = affected_clusters.rbegin(); it != affected_clusters.rend(); ++it) {
        Cluster::ptr clust = *it;
        
        // Get site-aligned position for this cluster
        double site_aligned_x = get_site_x(clust->x, sub_row_start_x, site_width);
        double x = site_aligned_x;
        
        for (int cell_idx : clust->member) {
            double dx = x - cells[cell_idx].original_x;
            double dy = sub_row_y - cells[cell_idx].original_y;
            double displacement = std::sqrt(dx * dx + dy * dy);
            
            if (displacement > max_displacement_constraint) {
                return false;
            }
            
            x += cells[cell_idx].width;
        }
    }
    
    return true;
}

PlacementResult PlaceRow::place_row_incremental_trial(SubRow* sub_row, int new_cell_index,
                                                     const std::vector<Cell>& cells, 
                                                     double max_displacement_constraint,
                                                     bool add_penalty) {
    PlacementResult result;
    result.valid = false;
    result.cost = std::numeric_limits<double>::infinity();
    
    if (!sub_row) {
        return result;
    }
    
    const Cell& new_cell = cells[new_cell_index];
    
    // Check free width first
    if (new_cell.width > sub_row->free_width) {
        return result;
    }
    
    // Adjust cell position to be within sub-row bounds
    double cell_x = new_cell.original_x;
    if (cell_x < sub_row->start_x) {
        cell_x = sub_row->start_x;
    } else if (cell_x > sub_row->end_x - new_cell.width) {
        cell_x = sub_row->end_x - new_cell.width;
    }
    
    // Work with a copy of the cluster chain for trial
    Cluster::ptr trial_last_cluster = nullptr;
    if (sub_row->last_cluster) {
        // Deep copy the cluster chain
        std::vector<Cluster::ptr> clusters;
        Cluster::ptr current = sub_row->last_cluster;
        
        while (current) {
            clusters.push_back(current);
            current = current->predecessor;
        }
        
        Cluster::ptr prev_cloned = nullptr;
        for (auto it = clusters.rbegin(); it != clusters.rend(); ++it) {
            auto cloned = std::make_shared<Cluster>((*it)->x, prev_cloned);
            cloned->total_weight = (*it)->total_weight;
            cloned->total_width = (*it)->total_width;
            cloned->q_value = (*it)->q_value;
            cloned->member = (*it)->member;
            prev_cloned = cloned;
        }
        trial_last_cluster = prev_cloned;
    }
    
    // Check if new cell overlaps with last cluster
    if (!trial_last_cluster || trial_last_cluster->x + trial_last_cluster->total_width <= cell_x) {
        // Create new independent cluster
        auto new_cluster = std::make_shared<Cluster>(cell_x, trial_last_cluster);
        new_cluster->member.push_back(new_cell_index);
        new_cluster->total_weight = new_cell.weight;
        new_cluster->q_value = new_cell.weight * cell_x;
        new_cluster->total_width = new_cell.width;
        new_cluster->x = cell_x;
        
        result.modified_last_cluster = new_cluster;
    } else {
        // Add cell to last cluster
        trial_last_cluster->member.push_back(new_cell_index);
        trial_last_cluster->total_weight += new_cell.weight;
        trial_last_cluster->q_value += new_cell.weight * (cell_x - trial_last_cluster->total_width);
        trial_last_cluster->total_width += new_cell.width;
        
        // Collapse clusters
        trial_last_cluster = collapse_cluster_incremental(trial_last_cluster, 
                                                         sub_row->start_x, 
                                                         sub_row->end_x,
                                                         sub_row->site_width);
        
        // Check constraints with site alignment if add_penalty is true
        if (add_penalty && !check_cluster_constraints_with_site_alignment(
                trial_last_cluster, cells, max_displacement_constraint, 
                sub_row->y, sub_row->start_x, sub_row->site_width, true)) {
            result.valid = false;
            return result;
        }
        
        result.modified_last_cluster = trial_last_cluster;
    }
    
    // Calculate cell positions with site alignment
    positions_from_clusters_with_site_alignment(result.modified_last_cluster, 
                                               result.cell_positions, 
                                               cells, 
                                               sub_row->start_x, 
                                               sub_row->site_width);
    
    // Validate result and calculate cost for new cell
    result.valid = true;
    for (const auto& pos : result.cell_positions) {
        int cell_idx = pos.first;
        double x = pos.second;
        
        // Check boundaries
        if (x < sub_row->start_x - EPSILON || 
            x + cells[cell_idx].width > sub_row->end_x + EPSILON) {
            result.valid = false;
            return result;
        }
        
        // Check displacement constraint
        double dx = x - cells[cell_idx].original_x;
        double dy = sub_row->y - cells[cell_idx].original_y;
        double displacement = std::sqrt(dx * dx + dy * dy);
        
        if (displacement > max_displacement_constraint) {
            result.valid = false;
            return result;
        }
        
        // Cost for new cell
        if (cell_idx == new_cell_index) {
            result.cost = displacement;
        }
    }
    
    return result;
}

void PlaceRow::place_row_incremental_final(SubRow* sub_row, int new_cell_index,
                                          std::vector<Cell>& cells,
                                          const PlacementResult& result) {
    if (!result.valid || !sub_row) {
        return;
    }
    
    const Cell& new_cell = cells[new_cell_index];
    
    // Update free width
    sub_row->free_width -= new_cell.width;
    
    // Adjust cell position to be within sub-row bounds
    double cell_x = new_cell.original_x;
    if (cell_x < sub_row->start_x) {
        cell_x = sub_row->start_x;
    } else if (cell_x > sub_row->end_x - new_cell.width) {
        cell_x = sub_row->end_x - new_cell.width;
    }
    
    // Update the actual cluster state
    if (!sub_row->last_cluster || sub_row->last_cluster->x + sub_row->last_cluster->total_width <= cell_x) {
        // Create new independent cluster
        sub_row->last_cluster = std::make_shared<Cluster>(cell_x, sub_row->last_cluster);
        sub_row->last_cluster->member.push_back(new_cell_index);
        sub_row->last_cluster->total_weight = new_cell.weight;
        sub_row->last_cluster->q_value = new_cell.weight * cell_x;
        sub_row->last_cluster->total_width = new_cell.width;
    } else {
        // Add cell to last cluster
        sub_row->last_cluster->member.push_back(new_cell_index);
        sub_row->last_cluster->total_weight += new_cell.weight;
        sub_row->last_cluster->q_value += new_cell.weight * (cell_x - sub_row->last_cluster->total_width);
        sub_row->last_cluster->total_width += new_cell.width;
        
        // Collapse clusters
        sub_row->last_cluster = collapse_cluster_incremental(sub_row->last_cluster,
                                                            sub_row->start_x,
                                                            sub_row->end_x,
                                                            sub_row->site_width);
    }
    
    // Apply positions to cells with site alignment
    for (const auto& pos : result.cell_positions) {
        cells[pos.first].current_x = pos.second;
        cells[pos.first].current_y = sub_row->y;
    }
}

void PlaceRow::positions_from_clusters_with_site_alignment(
    Cluster::ptr last_cluster,
    std::vector<std::pair<int, double>>& positions,
    const std::vector<Cell>& cells,
    double sub_row_start_x,
    int site_width) {
    
    positions.clear();
    
    // Collect all clusters
    std::vector<Cluster::ptr> clusters;
    Cluster::ptr current = last_cluster;
    while (current) {
        clusters.push_back(current);
        current = current->predecessor;
    }
    
    // Process clusters in forward order (reverse of collection order)
    for (auto it = clusters.rbegin(); it != clusters.rend(); ++it) {
        Cluster::ptr cluster = *it;
        
        // Get site-aligned position for cluster start
        double x = get_site_x(cluster->x, sub_row_start_x, site_width);
        
        for (int cell_idx : cluster->member) {
            positions.push_back({cell_idx, x});
            x += cells[cell_idx].width;
        }
    }
}

double PlaceRow::get_site_x(double x, double min_x, int site_width) {
    double shift_x = x - min_x;
    return min_x + std::round(shift_x / site_width) * site_width;
}