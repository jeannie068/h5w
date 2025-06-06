// Simplified place_row.cpp - matching reference implementation
#include "place_row.hpp"
#include <algorithm>
#include <cmath>
#include <stack>
#include <limits>

std::pair<int, double> PlaceRow::place_row_trial(SubRow* sub_row, 
                                                 Cell* new_cell,
                                                 const std::vector<Cell>& cells,
                                                 double max_displacement_constraint,
                                                 bool add_penalty) {
    // Check if cell can fit
    if (!sub_row || new_cell->width > sub_row->free_width) {
        return {-1, std::numeric_limits<double>::max()};
    }
    
    // Adjust cell position to be within sub-row bounds
    double cell_x = new_cell->original_x;
    if (cell_x < sub_row->start_x) {
        cell_x = sub_row->start_x;
    } else if (cell_x > sub_row->end_x - new_cell->width) {
        cell_x = sub_row->end_x - new_cell->width;
    }
    
    Cluster::ptr cluster = sub_row->last_cluster;
    
    // Simple case: no overlap with existing clusters
    if (!cluster || cluster->x + cluster->total_width <= cell_x) {
        new_cell->current_x = cell_x;
        new_cell->current_y = sub_row->y;
        
        if (add_penalty) {
            double site_x = get_site_x(cell_x, sub_row->start_x, sub_row->site_width);
            double dx = site_x - new_cell->original_x;
            double dy = sub_row->y - new_cell->original_y;
            double displacement = std::sqrt(dx * dx + dy * dy);
            
            if (displacement > max_displacement_constraint) {
                return {-1, std::numeric_limits<double>::max()};
            }
        }
        
        double dx = new_cell->current_x - new_cell->original_x;
        double dy = new_cell->current_y - new_cell->original_y;
        return {0, std::sqrt(dx * dx + dy * dy)};
    }
    
    // Complex case: need to merge with clusters
    // Add cell to cluster temporarily
    int cluster_weight = cluster->total_weight + new_cell->weight;
    double cluster_q = cluster->q_value + new_cell->weight * (cell_x - cluster->total_width);
    int cluster_width = cluster->total_width + new_cell->width;
    
    std::stack<Cluster::ptr> cluster_stack;
    double cluster_x = 0;
    
    // Collapse clusters
    while (true) {
        cluster_stack.push(cluster);
        cluster_x = cluster_q / cluster_weight;
        
        if (cluster_x < sub_row->start_x) {
            cluster_x = sub_row->start_x;
        }
        if (cluster_x > sub_row->end_x - cluster_width) {
            cluster_x = sub_row->end_x - cluster_width;
        }
        
        Cluster::ptr prev_cluster = cluster->predecessor;
        if (prev_cluster && prev_cluster->x + prev_cluster->total_width > cluster_x) {
            // Merge with predecessor
            cluster_q = prev_cluster->q_value + cluster_q - cluster_weight * prev_cluster->total_width;
            cluster_weight = prev_cluster->total_weight + cluster_weight;
            cluster_width = prev_cluster->total_width + cluster_width;
            cluster = prev_cluster;
        } else {
            break;
        }
    }
    
    // Calculate new cell position
    new_cell->current_x = cluster_x + cluster_width - new_cell->width;
    new_cell->current_y = sub_row->y;
    
    // Check constraints if add_penalty is true
    if (add_penalty) {
        // Check all cells in affected clusters
        int x = get_site_x(cluster_x, sub_row->start_x, sub_row->site_width);
        
        while (!cluster_stack.empty()) {
            for (int cell_idx : cluster_stack.top()->member) {
                const Cell& cell = cells[cell_idx];
                double dx = x - cell.original_x;
                double dy = sub_row->y - cell.original_y;
                if (std::sqrt(dx * dx + dy * dy) > max_displacement_constraint) {
                    return {-1, std::numeric_limits<double>::max()};
                }
                x += cell.width;
            }
            cluster_stack.pop();
        }
        
        // Check new cell
        double dx = x - new_cell->original_x;
        double dy = sub_row->y - new_cell->original_y;
        if (std::sqrt(dx * dx + dy * dy) > max_displacement_constraint) {
            return {-1, std::numeric_limits<double>::max()};
        }
    }
    
    // Return displacement
    double dx = new_cell->current_x - new_cell->original_x;
    double dy = new_cell->current_y - new_cell->original_y;
    return {0, std::sqrt(dx * dx + dy * dy)};
}

void PlaceRow::place_row_final(SubRow* sub_row, 
                              Cell* new_cell,
                              std::vector<Cell>& cells) {
    // Update free width
    sub_row->free_width -= new_cell->width;
    
    // Get cell index properly
    int cell_index = new_cell - &cells[0];
    
    // Adjust cell position to be within sub-row bounds
    double cell_x = new_cell->original_x;
    if (cell_x < sub_row->start_x) {
        cell_x = sub_row->start_x;
    } else if (cell_x > sub_row->end_x - new_cell->width) {
        cell_x = sub_row->end_x - new_cell->width;
    }
    
    Cluster::ptr cluster = sub_row->last_cluster;
    
    // Simple case: create new cluster
    if (!cluster || cluster->x + cluster->total_width <= cell_x) {
        sub_row->last_cluster = std::make_shared<Cluster>(cell_x, cluster);
        cluster = sub_row->last_cluster;
        
        cluster->member.push_back(cell_index);
        cluster->total_weight = new_cell->weight;
        cluster->q_value = new_cell->weight * cell_x;
        cluster->total_width = new_cell->width;
    } else {
        // Add cell to existing cluster
        cluster->member.push_back(cell_index);
        cluster->total_weight += new_cell->weight;
        cluster->q_value += new_cell->weight * (cell_x - cluster->total_width);
        cluster->total_width += new_cell->width;
        
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
}

double PlaceRow::get_site_x(double x, double min_x, int site_width) {
    double shift_x = x - min_x;
    return min_x + std::round(shift_x / site_width) * site_width;
}