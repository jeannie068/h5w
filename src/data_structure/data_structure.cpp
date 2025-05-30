#include "data_structure.hpp"
#include <algorithm>
#include <cmath>

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
    
    // Check if there's vertical overlap
    return !(blockage_top <= row_y || y >= row_top);
}

bool Blockage::overlaps_row_horizontally(int start_x, int end_x) const {
    int blockage_right = x + width;
    
    // Check if there's horizontal overlap
    return !(blockage_right <= start_x || x >= end_x);
}

// SubRow methods
bool SubRow::can_fit_cell(int cell_width) const {
    int required_sites = static_cast<int>(std::ceil(static_cast<double>(cell_width) / site_width));
    int current_used_sites = 0;
    
    // Calculate currently used sites (this is approximate, actual calculation needs cell positions)
    return required_sites <= available_sites;
}

double SubRow::get_site_aligned_x(double x) const {
    // Align x to site boundary using floor (align to left or equal site boundary)
    double relative_x = x - start_x;
    int site_offset = static_cast<int>(std::floor(relative_x / site_width));
    double aligned_x = start_x + site_offset * site_width;
    
    // Ensure it's not less than sub-row start due to floating point issues
    aligned_x = std::max(start_x, aligned_x);
    
    return aligned_x;
}

void SubRow::add_cell(int cell_index, const std::vector<Cell>& all_cells) {
    // Insert cell maintaining sorted order by original_x
    auto insert_pos = std::lower_bound(cells.begin(), cells.end(), cell_index,
        [&all_cells](int a, int b) {
            return all_cells[a].original_x < all_cells[b].original_x;
        });
    cells.insert(insert_pos, cell_index);
}

void SubRow::remove_cell(int cell_index) {
    auto it = std::find(cells.begin(), cells.end(), cell_index);
    if (it != cells.end()) {
        cells.erase(it);
    }
}

// Row methods
void Row::create_sub_rows(const std::vector<Blockage>& blockages, int row_index) {
    std::vector<std::pair<int, int>> blockage_intervals;
    
    // Find all blockages that overlap with this row
    for (const auto& blockage : blockages) {
        if (blockage.overlaps_row_vertically(y, height) && 
            blockage.overlaps_row_horizontally(start_x, end_x)) {
            int block_start = std::max(blockage.x, start_x);
            int block_end = std::min(blockage.x + blockage.width, end_x);
            blockage_intervals.push_back({block_start, block_end});
        }
    }
    
    // Sort blockage intervals by start position
    std::sort(blockage_intervals.begin(), blockage_intervals.end());
    
    // Merge overlapping intervals
    std::vector<std::pair<int, int>> merged_intervals;
    for (const auto& interval : blockage_intervals) {
        if (merged_intervals.empty() || merged_intervals.back().second < interval.first) {
            merged_intervals.push_back(interval);
        } else {
            merged_intervals.back().second = std::max(merged_intervals.back().second, interval.second);
        }
    }
    
    // Create sub-rows in the gaps
    sub_rows.clear();
    int current_x = start_x;
    
    for (const auto& interval : merged_intervals) {
        if (current_x < interval.first) {
            // Create sub-row before this blockage
            double sub_start = current_x;
            double sub_end = interval.first;
            
            // Ensure sub-row boundaries are site-aligned
            // sub_start: ceil to next site boundary
            double relative_start = sub_start - start_x;
            int start_site_offset = static_cast<int>(std::ceil(relative_start / site_width));
            sub_start = start_x + start_site_offset * site_width;
            
            // sub_end: floor to previous site boundary  
            double relative_end = sub_end - start_x;
            int end_site_offset = static_cast<int>(std::floor(relative_end / site_width));
            sub_end = start_x + end_site_offset * site_width;
            
            if (sub_end > sub_start) {
                sub_rows.emplace_back(row_index, sub_start, sub_end, y, height, site_width);
            }
        }
        current_x = interval.second;
    }
    
    // Create sub-row after the last blockage (if any space left)
    if (current_x < end_x) {
        double sub_start = current_x;
        double sub_end = end_x;
        
        // Ensure sub-row boundaries are site-aligned
        // sub_start: ceil to next site boundary
        double relative_start = sub_start - start_x;
        int start_site_offset = static_cast<int>(std::ceil(relative_start / site_width));
        sub_start = start_x + start_site_offset * site_width;
        
        // sub_end: should already be site-aligned (end_x = start_x + site_width * num_sites)
        // but apply floor for consistency
        double relative_end = sub_end - start_x;
        int end_site_offset = static_cast<int>(std::floor(relative_end / site_width));
        sub_end = start_x + end_site_offset * site_width;
        
        if (sub_end > sub_start) {
            sub_rows.emplace_back(row_index, sub_start, sub_end, y, height, site_width);
        }
    }
    
    // If no blockages overlap, create one sub-row for the entire row
    if (merged_intervals.empty()) {
        sub_rows.emplace_back(row_index, start_x, end_x, y, height, site_width);
    }
}

// PlacementData methods
void PlacementData::initialize_sub_row_pointers() {
    all_sub_rows.clear();
    for (auto& row : rows) {
        for (auto& sub_row : row.sub_rows) {
            all_sub_rows.push_back(&sub_row);
        }
    }
}

bool PlacementData::violates_max_displacement(int cell_index) const {
    if (cell_index < 0 || cell_index >= static_cast<int>(cells.size())) {
        return true;
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
    for (const auto& cell : cells) {
        max_disp = std::max(max_disp, cell.get_displacement());
    }
    return std::ceil(max_disp);
}

int PlacementData::get_sub_row_containing_cell(int cell_index) const {
    for (int i = 0; i < static_cast<int>(all_sub_rows.size()); ++i) {
        const auto& sub_row = *all_sub_rows[i];
        if (std::find(sub_row.cells.begin(), sub_row.cells.end(), cell_index) != sub_row.cells.end()) {
            return i;
        }
    }
    return -1;  // Cell not found in any sub-row
}