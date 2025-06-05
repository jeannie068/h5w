// data_structure.hpp - Modified version with incremental cluster support
#ifndef DATA_STRUCTURE_HPP
#define DATA_STRUCTURE_HPP

#include <vector>
#include <string>
#include <cmath>
#include <memory>

struct Cell {
    std::string name;
    int width;
    int height;
    double original_x;
    double original_y;
    double current_x;
    double current_y;
    double weight;  // e(i) in paper, default to 1.0
    bool is_legalized;
    
    Cell(const std::string& n, int w, int h, double ox, double oy)
        : name(n), width(w), height(h), original_x(ox), original_y(oy),
          current_x(ox), current_y(oy), weight(1.0), is_legalized(false) {}
    
    // Calculate Euclidean displacement
    double get_displacement() const;
};

// Forward declaration
struct Cluster;

struct Cluster {
    using ptr = std::shared_ptr<Cluster>;
    
    double x;                    // x position of cluster
    double total_weight;         // e_c
    double total_width;          // w_c
    double q_value;              // q_c
    ptr predecessor;             // Link to previous cluster
    std::vector<int> member;     // Cell indices in this cluster
    
    Cluster() : x(0.0), total_weight(0.0), total_width(0.0), q_value(0.0) {}
    
    Cluster(double x_pos, ptr pred) 
        : x(x_pos), total_weight(0.0), total_width(0.0), q_value(0.0), predecessor(pred) {}
};

struct Blockage {
    std::string name;
    int width;
    int height;
    int x;
    int y;
    
    Blockage(const std::string& n, int w, int h, int bx, int by)
        : name(n), width(w), height(h), x(bx), y(by) {}
    
    // Check if blockage overlaps with row vertically
    bool overlaps_row_vertically(int row_y, int row_height) const;
    
    // Check if blockage overlaps with row horizontally in given range
    bool overlaps_row_horizontally(int start_x, int end_x) const;
};

struct SubRow {
    int parent_row_index;
    double start_x;
    double end_x;
    int y;
    int height;
    int site_width;
    std::vector<int> cells;  // indices of cells in this sub-row
    int available_sites;
    
    // For incremental cluster management
    Cluster::ptr last_cluster;   // The rightmost cluster in this sub-row
    
    SubRow(int parent_idx, double sx, double ex, int row_y, int row_h, int sw)
        : parent_row_index(parent_idx), start_x(sx), end_x(ex), 
          y(row_y), height(row_h), site_width(sw), last_cluster(nullptr) {
        available_sites = static_cast<int>((end_x - start_x) / site_width);
    }
    
    // Check if a cell can fit in this sub-row
    bool can_fit_cell(int cell_width) const;
    
    // Get site-aligned x position
    double get_site_aligned_x(double x) const;
    
    // Add cell to this sub-row (maintains sorted order by original_x)
    void add_cell(int cell_index, const std::vector<Cell>& all_cells);
    
    // Remove cell from this sub-row
    void remove_cell(int cell_index);
    
    // Clone this sub-row for trial purposes
    std::unique_ptr<SubRow> clone() const;
};

struct Row {
    std::string name;
    int site_width;
    int height;
    int row_start_x;
    int y;
    int num_sites;
    int row_end_x;  // calculated: start_x + site_width * num_sites
    std::vector<SubRow> sub_rows;
    
    Row(const std::string& n, int sw, int h, int sx, int row_y, int ns)
        : name(n), site_width(sw), height(h), row_start_x(sx), y(row_y), num_sites(ns) {
        row_end_x = row_start_x + site_width * num_sites;
    }
    
    // Create sub-rows by splitting around blockages
    void create_sub_rows(const std::vector<Blockage>& blockages, int row_index);
};

struct PlacementData {
    double max_displacement_constraint;
    std::vector<Cell> cells;
    std::vector<Blockage> blockages;
    std::vector<Row> rows;
    std::vector<SubRow*> all_sub_rows;  // pointers to all sub-rows for easy access
    
    PlacementData() : max_displacement_constraint(0.0) {}
    
    // Initialize all_sub_rows pointers after creating sub-rows
    void initialize_sub_row_pointers();
    
    // Check if a cell violates max displacement constraint
    bool violates_max_displacement(int cell_index) const;
    
    // Calculate total displacement (sum of all cell displacements)
    double calculate_total_displacement() const;
    
    // Calculate max displacement (maximum among all cells)
    double calculate_max_displacement() const;
    
    // Get the sub-row index that contains the given cell
    int get_sub_row_containing_cell(int cell_index) const;
};

#endif // DATA_STRUCTURE_HPP