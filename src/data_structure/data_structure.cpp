#include "data_structure.hpp"
// #include "../Logger.hpp" // Assuming Logger might be used, uncomment if needed
#include <algorithm>
#include <cmath>
#include <iostream> // For potential debugging, otherwise can be removed

// Cell methods
double Cell::get_displacement() const {
    double dx = current_x - original_x;
    double dy = current_y - original_y;
    return std::sqrt(dx * dx + dy * dy);
}

// Blockage methods
bool Blockage::overlaps_row_vertically(int row_y, int row_height) const {
    int blockage_top = y + height;
    int row_top = row_y + row_height;
    return !(blockage_top <= row_y || y >= row_top);
}

bool Blockage::overlaps_row_horizontally(int start_x, int end_x) const {
    int blockage_right = x + width;
    return !(blockage_right <= start_x || x >= end_x);
}

// SubRow methods
bool SubRow::can_fit_cell(int cell_width) const {
    // This is a simplified check. A more accurate check would consider remaining_sites.
    // For now, ensure the sub-row itself is wide enough.
    // The available_sites check is done in AbacusLegalizer::try_place_cell_trial
    return (end_x - start_x) >= cell_width;
}

double SubRow::get_site_aligned_x(double x_coord) const {
    if (site_width == 0) return x_coord; // Avoid division by zero
    double relative_x = x_coord - start_x;
    // Ensure alignment doesn't push cell before subrow_start_x due to floating point
    double aligned_x = start_x + std::round(relative_x / site_width) * site_width;
    return std::max(start_x, aligned_x);
}

void SubRow::add_cell_to_final_list(int cell_index, const std::vector<Cell>& all_cells) {
    auto it = std::lower_bound(cells.begin(), cells.end(), cell_index, 
        [&](int a_idx, int b_val_idx) {
            return all_cells[a_idx].original_x < all_cells[b_val_idx].original_x;
        });
    cells.insert(it, cell_index);
}

void SubRow::remove_cell_from_final_list(int cell_index) {
    auto it = std::find(cells.begin(), cells.end(), cell_index);
    if (it != cells.end()) {
        cells.erase(it);
    }
}


// Row methods
void Row::create_sub_rows(const std::vector<Blockage>& blockages_vec, int row_index) {
    std::vector<std::pair<int, int>> blockage_intervals;
    for (const auto& blockage : blockages_vec) {
        if (blockage.overlaps_row_vertically(y, height) &&
            blockage.overlaps_row_horizontally(row_start_x, row_end_x)) {
            blockage_intervals.push_back({std::max(blockage.x, row_start_x),
                                          std::min(blockage.x + blockage.width, row_end_x)});
        }
    }

    std::sort(blockage_intervals.begin(), blockage_intervals.end());

    std::vector<std::pair<int, int>> merged_intervals;
    for (const auto& interval : blockage_intervals) {
        if (merged_intervals.empty() || merged_intervals.back().second < interval.first) {
            merged_intervals.push_back(interval);
        } else {
            merged_intervals.back().second = std::max(merged_intervals.back().second, interval.second);
        }
    }

    sub_rows.clear();
    double current_s_x = static_cast<double>(row_start_x);

    for (const auto& interval : merged_intervals) {
        if (current_s_x < static_cast<double>(interval.first)) {
            double sub_start = current_s_x;
            double sub_end = static_cast<double>(interval.first);
             // Ensure sub-row boundaries are site-aligned for start
            double rel_start = sub_start - row_start_x;
            int start_site_offset = static_cast<int>(std::ceil(rel_start / site_width));
            sub_start = static_cast<double>(row_start_x) + start_site_offset * site_width;

            // Ensure sub-row boundaries are site-aligned for end
            double rel_end = sub_end - row_start_x;
            int end_site_offset = static_cast<int>(std::floor(rel_end / site_width));
            sub_end = static_cast<double>(row_start_x) + end_site_offset * site_width;


            if (sub_end > sub_start + 1e-9) { // Add epsilon for float comparison
                 sub_rows.emplace_back(row_index, sub_start, sub_end, y, height, site_width);
            }
        }
        current_s_x = static_cast<double>(interval.second);
    }

    if (current_s_x < static_cast<double>(row_end_x)) {
        double sub_start = current_s_x;
        double sub_end = static_cast<double>(row_end_x);

        double rel_start = sub_start - row_start_x;
        int start_site_offset = static_cast<int>(std::ceil(rel_start / site_width));
        sub_start = static_cast<double>(row_start_x) + start_site_offset * site_width;
        
        // row_end_x is already site aligned: row_start_x + num_sites * site_width
        // No need to realign sub_end if it's row_end_x

        if (sub_end > sub_start + 1e-9) {
            sub_rows.emplace_back(row_index, sub_start, sub_end, y, height, site_width);
        }
    }
}

// PlacementData methods
void PlacementData::initialize_sub_row_pointers() {
    all_sub_rows.clear();
    for (size_t i = 0; i < rows.size(); ++i) {
        for (size_t j = 0; j < rows[i].sub_rows.size(); ++j) {
            all_sub_rows.push_back(&rows[i].sub_rows[j]);
        }
    }
     // Assign IDs to cells
    for(size_t i=0; i < cells.size(); ++i) {
        cells[i].id = i;
    }
}

bool PlacementData::violates_max_displacement(int cell_index) const {
    if (cell_index < 0 || cell_index >= static_cast<int>(cells.size())) {
        return true; // Or handle error
    }
    return cells[cell_index].get_displacement() > max_displacement_constraint;
}

double PlacementData::calculate_total_displacement() const {
    double total = 0.0;
    for (const auto& cell : cells) {
        total += cell.get_displacement();
    }
    return std::ceil(total);
}

double PlacementData::calculate_max_displacement() const {
    double max_disp = 0.0;
    if (cells.empty()) return 0.0;
    for (const auto& cell : cells) {
        max_disp = std::max(max_disp, cell.get_displacement());
    }
    return std::ceil(max_disp);
}

int PlacementData::get_sub_row_containing_cell(int cell_index) const {
    for (size_t i = 0; i < all_sub_rows.size(); ++i) {
        const auto& sub_row_cells = all_sub_rows[i]->cells;
        if (std::find(sub_row_cells.begin(), sub_row_cells.end(), cell_index) != sub_row_cells.end()) {
            return i;
        }
    }
    return -1;
}
