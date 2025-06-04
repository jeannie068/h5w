#ifndef DATA_STRUCTURE_HPP
#define DATA_STRUCTURE_HPP

#include <vector>
#include <string>
#include <cmath>
#include <memory> // Required for std::shared_ptr
#include <algorithm> // Required for std::sort, std::find

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
    int id; // Unique ID for cells, can be its index in the main cells vector

    Cell(const std::string& n, int w, int h, double ox, double oy, int cell_id = -1)
        : name(n), width(w), height(h), original_x(ox), original_y(oy),
          current_x(ox), current_y(oy), weight(1.0), is_legalized(false), id(cell_id) {}

    double get_displacement() const;
};

struct Cluster {
    double total_weight;    // e_c
    double total_width;     // w_c
    double q_value;         // q_c
    double optimal_x;       // calculated optimal x position of the cluster's left edge
    
    std::vector<int> member_cell_indices; // Indices of cells in this cluster, sorted by original_x
    std::shared_ptr<Cluster> predecessor;

    Cluster() : total_weight(0.0), total_width(0.0), q_value(0.0), optimal_x(0.0), predecessor(nullptr) {}

    // Constructor for a new cluster with a single cell
    Cluster(int cell_idx, const Cell& cell_data, double initial_x, std::shared_ptr<Cluster> pred)
        : total_weight(cell_data.weight),
          total_width(cell_data.width),
          q_value(cell_data.weight * initial_x), // q_c = e(i) * (x'(i) - 0) since it's the first cell in this cluster context
          optimal_x(initial_x),
          predecessor(pred) {
        member_cell_indices.push_back(cell_idx);
    }
};

struct Blockage {
    std::string name;
    int width;
    int height;
    int x;
    int y;

    Blockage(const std::string& n, int w, int h, int bx, int by)
        : name(n), width(w), height(h), x(bx), y(by) {}

    bool overlaps_row_vertically(int row_y, int row_height) const;
    bool overlaps_row_horizontally(int start_x, int end_x) const;
};

struct SubRow {
    int parent_row_index;
    double start_x;
    double end_x;
    int y;
    int height;
    int site_width;
    std::vector<int> cells;  // Finalized indices of cells in this sub-row, kept sorted by original_x
    int available_sites;

    std::shared_ptr<Cluster> last_cluster; // Points to the rightmost cluster in this sub-row

    SubRow(int parent_idx, double sx, double ex, int row_y, int row_h, int sw)
        : parent_row_index(parent_idx), start_x(sx), end_x(ex),
          y(row_y), height(row_h), site_width(sw), last_cluster(nullptr) {
        available_sites = static_cast<int>((end_x - start_x) / site_width);
    }

    bool can_fit_cell(int cell_width) const;
    double get_site_aligned_x(double x_coord) const;
    
    void add_cell_to_final_list(int cell_index, const std::vector<Cell>& all_cells);
    void remove_cell_from_final_list(int cell_index);
};

struct Row {
    std::string name;
    int site_width;
    int height;
    int row_start_x;
    int y;
    int num_sites;
    int row_end_x;
    std::vector<SubRow> sub_rows;

    Row(const std::string& n, int sw, int h, int sx, int row_y, int ns)
        : name(n), site_width(sw), height(h), row_start_x(sx), y(row_y), num_sites(ns) {
        row_end_x = row_start_x + site_width * num_sites;
    }

    void create_sub_rows(const std::vector<Blockage>& blockages, int row_index);
};

struct PlacementData {
    double max_displacement_constraint;
    std::vector<Cell> cells;
    std::vector<Blockage> blockages;
    std::vector<Row> rows;
    std::vector<SubRow*> all_sub_rows;

    PlacementData() : max_displacement_constraint(0.0) {}

    void initialize_sub_row_pointers();
    bool violates_max_displacement(int cell_index) const;
    double calculate_total_displacement() const;
    double calculate_max_displacement() const;
    int get_sub_row_containing_cell(int cell_index) const;
};

#endif // DATA_STRUCTURE_HPP
