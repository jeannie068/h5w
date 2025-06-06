// place_row.cpp - Modified key functions
#include "place_row.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stack>

// Fixed place_row_incremental_trial in place_row.cpp

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
    
    // Quick check: free width first
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
    
    // Simple case: no existing clusters or new cell doesn't overlap
    if (!sub_row->last_cluster || sub_row->last_cluster->x + sub_row->last_cluster->total_width <= cell_x) {
        // Cell will be placed at cell_x (after site alignment)
        double site_aligned_x = get_site_x(cell_x, sub_row->start_x, sub_row->site_width);
        
        // Calculate displacement
        double dx = site_aligned_x - new_cell.original_x;
        double dy = sub_row->y - new_cell.original_y;
        double displacement = std::sqrt(dx * dx + dy * dy);
        
        // Check constraint if add_penalty is true
        if (add_penalty && displacement > max_displacement_constraint) {
            return result;
        }
        
        result.cost = displacement;
        result.valid = true;
        result.cell_positions.push_back({new_cell_index, site_aligned_x});
        return result;
    }
    
    // Complex case: need to handle cluster collapse
    // Create temporary cluster state (following reference implementation approach)
    std::vector<Cluster::ptr> cluster_stack;
    Cluster::ptr current = sub_row->last_cluster;
    
    // Simulate adding cell to last cluster
    int cluster_weight = current->total_weight + new_cell.weight;
    double cluster_q = current->q_value + new_cell.weight * (cell_x - current->total_width);
    int cluster_width = current->total_width + new_cell.width;
    
    // Collapse clusters (following reference implementation logic)
    while (true) {
        cluster_stack.push_back(current);
        
        double cluster_x = cluster_q / cluster_weight;
        if (cluster_x < sub_row->start_x) {
            cluster_x = sub_row->start_x;
        }
        if (cluster_x > sub_row->end_x - cluster_width) {
            cluster_x = sub_row->end_x - cluster_width;
        }
        
        Cluster::ptr prev = current->predecessor;
        if (prev && prev->x + prev->total_width > cluster_x) {
            // Merge with predecessor
            cluster_q = prev->q_value + cluster_q - cluster_weight * prev->total_width;
            cluster_weight = prev->total_weight + cluster_weight;
            cluster_width = prev->total_width + cluster_width;
            current = prev;
        } else {
            // No more collapsing needed
            // Now calculate positions and check constraints
            
            // Get all clusters that are not affected
            std::vector<Cluster::ptr> unaffected_clusters;
            Cluster::ptr temp = sub_row->last_cluster;
            while (temp && temp != current) {
                temp = temp->predecessor;
            }
            if (temp) {
                temp = temp->predecessor;
                while (temp) {
                    unaffected_clusters.push_back(temp);
                    temp = temp->predecessor;
                }
            }
            
            // Calculate positions for all cells
            result.cell_positions.clear();
            
            // First, unaffected clusters (in reverse order)
            for (auto it = unaffected_clusters.rbegin(); it != unaffected_clusters.rend(); ++it) {
                double x = get_site_x((*it)->x, sub_row->start_x, sub_row->site_width);
                for (int cell_idx : (*it)->member) {
                    result.cell_positions.push_back({cell_idx, x});
                    x += cells[cell_idx].width;
                }
            }
            
            // Then, affected clusters
            double x = get_site_x(cluster_x, sub_row->start_x, sub_row->site_width);
            
            // Process clusters from bottom to top of stack
            for (auto it = cluster_stack.rbegin(); it != cluster_stack.rend(); ++it) {
                for (int cell_idx : (*it)->member) {
                    result.cell_positions.push_back({cell_idx, x});
                    
                    if (add_penalty) {
                        double dx = x - cells[cell_idx].original_x;
                        double dy = sub_row->y - cells[cell_idx].original_y;
                        double disp = std::sqrt(dx * dx + dy * dy);
                        if (disp > max_displacement_constraint) {
                            // Constraint violated
                            result.valid = false;
                            return result;
                        }
                    }
                    
                    x += cells[cell_idx].width;
                }
            }
            
            // Add new cell
            result.cell_positions.push_back({new_cell_index, x});
            double dx = x - new_cell.original_x;
            double dy = sub_row->y - new_cell.original_y;
            double displacement = std::sqrt(dx * dx + dy * dy);
            
            if (add_penalty && displacement > max_displacement_constraint) {
                result.valid = false;
                return result;
            }
            
            result.cost = displacement;
            result.valid = true;
            break;
        }
    }
    
    return result;
}

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

// Fixed place_row_incremental_final in place_row.cpp

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
    
    // Update the actual cluster state (following reference's placeRowFinal)
    Cluster::ptr cluster = sub_row->last_cluster;
    
    if (!cluster || cluster->x + cluster->total_width <= cell_x) {
        // Create new independent cluster
        sub_row->last_cluster = std::make_shared<Cluster>(cell_x, cluster);
        cluster = sub_row->last_cluster;
        
        // Add cell
        cluster->member.push_back(new_cell_index);
        cluster->total_weight = new_cell.weight;
        cluster->q_value = new_cell.weight * cell_x;
        cluster->total_width = new_cell.width;
    } else {
        // Add cell to last cluster
        cluster->member.push_back(new_cell_index);
        cluster->total_weight += new_cell.weight;
        cluster->q_value += new_cell.weight * (cell_x - cluster->total_width);
        cluster->total_width += new_cell.width;
        
        // Collapse clusters
        while (true) {
            cluster->x = cluster->q_value / cluster->total_weight;
            
            if (cluster->x < sub_row->start_x) {
                cluster->x = sub_row->start_x;
            }
            if (cluster->x > sub_row->end_x - cluster->total_width) {
                cluster->x = sub_row->end_x - cluster->total_width;
            }
            
            Cluster::ptr prev_cluster = cluster->predecessor;
            if (prev_cluster && prev_cluster->x + prev_cluster->total_width > cluster->x) {
                // Merge clusters
                prev_cluster->member.insert(prev_cluster->member.end(), 
                                          cluster->member.begin(), 
                                          cluster->member.end());
                prev_cluster->total_weight += cluster->total_weight;
                prev_cluster->q_value += cluster->q_value - cluster->total_weight * prev_cluster->total_width;
                prev_cluster->total_width += cluster->total_width;
                
                cluster = prev_cluster;
            } else {
                break;
            }
        }
        
        sub_row->last_cluster = cluster;
    }
    
    // Note: We don't update cell positions here - that's done in determine_final_positions
}


double PlaceRow::get_site_x(double x, double min_x, int site_width) {
    // Following reference implementation's getSiteX
    double shift_x = x - min_x;
    return min_x + std::round(shift_x / site_width) * site_width;
}