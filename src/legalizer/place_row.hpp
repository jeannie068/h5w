#ifndef PLACE_ROW_HPP
#define PLACE_ROW_HPP

#include "../data_structure/data_structure.hpp"
#include <vector>

// Structure to hold placement result for trial mode
struct PlacementResult {
    bool valid;
    double cost;  // Displacement of the new cell
    std::vector<std::pair<int, double>> cell_positions;  // cell_index, x_position pairs
};

class PlaceRow {
public:
    // Floating point comparison epsilon
    static constexpr double EPSILON = 1e-9;
    
    // Trial mode - calculates positions without modifying cells
    static PlacementResult place_row_trial(SubRow* sub_row, int new_cell_index, 
                                          const std::vector<Cell>& cells, 
                                          double max_displacement_constraint);
    
    // Final mode - actually places cells based on calculated positions
    static void place_row_final(SubRow* sub_row, const PlacementResult& result,
                               std::vector<Cell>& cells);
    
    // Legacy function for backward compatibility
    static bool place_row(SubRow* sub_row, std::vector<Cell>& cells, 
                          double max_displacement_constraint, bool trial_mode = false);

private:
    // Helper function for trial mode cluster collapse
    static void collapse_clusters_trial(std::vector<Cluster>& clusters, int cluster_index,
                                       double sub_row_start_x, double sub_row_end_x);
    
    // Cluster operations (from original implementation)
    static void add_cell_to_cluster(Cluster& cluster, int cell_index, 
                                   const std::vector<Cell>& cells);
    
    static void add_cluster_to_cluster(Cluster& target_cluster, const Cluster& source_cluster);    
    
    static void sort_cells_by_original_x(std::vector<int>& cell_indices, 
                                       const std::vector<Cell>& cells);
};

#endif // PLACE_ROW_HPP