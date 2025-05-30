#ifndef PLACE_ROW_HPP
#define PLACE_ROW_HPP

#include "../data_structure/data_structure.hpp"
#include <vector>

class PlaceRow {
public:
    // Floating point comparison epsilon
    static constexpr double EPSILON = 1e-9;
    
    // Main PlaceRow function - implements Algorithm 2 from the paper
    // Returns true if placement is successful, false if constraints are violated
    static bool place_row(SubRow* sub_row, std::vector<Cell>& cells, 
                          double max_displacement_constraint, bool trial_mode = false);

private:
    // Cluster operations (based on Algorithm 2 functions)
    static void add_cell_to_cluster(Cluster& cluster, int cell_index, 
                                   const std::vector<Cell>& cells);
    
    static void add_cluster_to_cluster(Cluster& target_cluster, const Cluster& source_cluster);
    
    static bool collapse_cluster(std::vector<Cluster>& clusters, int cluster_index,
                               double sub_row_start_x, double sub_row_end_x, int site_width);
    
    // Calculate optimal position for a cluster (equation 5 from paper)
    static double calculate_cluster_optimal_x(const Cluster& cluster, 
                                             const std::vector<Cell>& cells);
    
    // Site alignment functions
    static double get_site_aligned_position(double x, double sub_row_start_x, int site_width);
    
    static bool apply_site_alignment(std::vector<Cluster>& clusters, SubRow* sub_row);
    
    // Constraint checking
    static bool check_max_displacement_constraint(const std::vector<int>& cell_indices,
                                                 const std::vector<Cell>& cells,
                                                 double max_displacement_constraint);
    
    static bool check_max_displacement_constraint_with_positions(const std::vector<int>& cell_indices,
                                                               const std::vector<Cell>& cells,
                                                               const std::vector<std::pair<double, double>>& positions,
                                                               double max_displacement_constraint);
    
    static bool check_row_boundary_constraint(const std::vector<Cluster>& clusters,
                                             double sub_row_start_x, double sub_row_end_x);
    
    // Helper functions
    static void transform_clusters_to_cell_positions(const std::vector<Cluster>& clusters,
                                                   SubRow* sub_row, std::vector<Cell>& cells);
    
    static double calculate_total_cluster_width(const Cluster& cluster);
    
    // Debug and validation functions
    static bool validate_non_overlap(const std::vector<int>& cell_indices,
                                   const std::vector<Cell>& cells);
    
    static void sort_cells_by_original_x(std::vector<int>& cell_indices, 
                                       const std::vector<Cell>& cells);
};

#endif // PLACE_ROW_HPP