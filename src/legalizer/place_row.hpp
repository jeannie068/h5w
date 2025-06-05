// place_row.hpp - Modified version
#ifndef PLACE_ROW_HPP
#define PLACE_ROW_HPP

#include "../data_structure/data_structure.hpp"
#include <vector>

// Structure to hold placement result for trial mode
struct PlacementResult {
    bool valid;
    double cost;  // Displacement of the new cell
    std::vector<std::pair<int, double>> cell_positions;  // cell_index, x_position pairs
    
    // For incremental clustering, store the modified cluster state
    Cluster::ptr modified_last_cluster;  
};

class PlaceRow {
public:
    // Floating point comparison epsilon
    static constexpr double EPSILON = 1e-9;
    
    // Incremental trial mode - calculates positions using existing clusters
    static PlacementResult place_row_incremental_trial(SubRow* sub_row, int new_cell_index, 
                                                      const std::vector<Cell>& cells, 
                                                      double max_displacement_constraint,
                                                      bool add_penalty = true);
    
    // Final mode - actually places cells and updates cluster state
    static void place_row_incremental_final(SubRow* sub_row, int new_cell_index,
                                          std::vector<Cell>& cells,
                                          const PlacementResult& result);
    
    // Get site-aligned x position
    static double get_site_x(double x, double min_x, int site_width);

private:
    // Helper function for incremental cluster collapse
    static Cluster::ptr collapse_cluster_incremental(Cluster::ptr cluster,
                                                    double sub_row_start_x, 
                                                    double sub_row_end_x,
                                                    int site_width);
    
    // Check if cluster violates constraints after collapse with site alignment
    static bool check_cluster_constraints_with_site_alignment(
        Cluster::ptr cluster,
        const std::vector<Cell>& cells,
        double max_displacement_constraint,
        int sub_row_y,
        double sub_row_start_x,
        int site_width,
        bool add_penalty);
    
    // Convert clusters to cell positions with site alignment
    static void positions_from_clusters_with_site_alignment(
        Cluster::ptr last_cluster,
        std::vector<std::pair<int, double>>& positions,
        const std::vector<Cell>& cells,
        double sub_row_start_x,
        int site_width);
};

#endif // PLACE_ROW_HPP