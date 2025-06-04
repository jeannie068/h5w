#include "data_structure.hpp"
#include "../Logger.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>

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
    // Logger::log("Creating sub-rows for row " + name + " (index=" + std::to_string(row_index) 
    //          + ", y=" + std::to_string(y) + ", x_range=[" + std::to_string(row_start_x) 
    //          + ", " + std::to_string(row_end_x) + "])");
    
    std::vector<std::pair<int, int>> blockage_intervals;
    
    // Find all blockages that overlap with this row
    for (const auto& blockage : blockages) {
        if (blockage.overlaps_row_vertically(y, height) && 
            blockage.overlaps_row_horizontally(row_start_x, row_end_x)) {
            int block_start = std::max(blockage.x, row_start_x);
            int block_end = std::min(blockage.x + blockage.width, row_end_x);
            blockage_intervals.push_back({block_start, block_end});
            
            // Logger::log("  Blockage " + blockage.name + " overlaps: [" 
            //          + std::to_string(block_start) + ", " + std::to_string(block_end) + "]");
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
    
    // std::cout << "  Merged blockage intervals: ";
    // Logger::log("  Merged blockage intervals: ");
    for (const auto& interval : merged_intervals) {
        // Logger::log("[" + std::to_string(interval.first) + ", " + std::to_string(interval.second) + "] ");
    }
    // std::cout << std::endl;
    
    // Create sub-rows in the gaps
    sub_rows.clear();
    int current_x = row_start_x;
    
    for (const auto& interval : merged_intervals) {
        if (current_x < interval.first) {
            // Create sub-row before this blockage
            double sub_start = current_x;
            double sub_end = interval.first;
            
            // Logger::log("  Trying to create sub-row before blockage: [" 
            //          + std::to_string(sub_start) + ", " + std::to_string(sub_end) + "]");
            
            // Ensure sub-row boundaries are site-aligned
            // sub_start: ceil to next site boundary
            double relative_start = sub_start - row_start_x;
            int start_site_offset = static_cast<int>(std::ceil(relative_start / site_width));
            sub_start = row_start_x + start_site_offset * site_width;
            
            // sub_end: floor to previous site boundary  
            double relative_end = sub_end - row_start_x;
            int end_site_offset = static_cast<int>(std::floor(relative_end / site_width));
            sub_end = row_start_x + end_site_offset * site_width;
            
            // Logger::log("    After site alignment: [" 
                    //  + std::to_string(sub_start) + ", " + std::to_string(sub_end) + "]");
            
            if (sub_end > sub_start) {
                sub_rows.emplace_back(row_index, sub_start, sub_end, y, height, site_width);
                // Logger::log(" - Created sub-row " + std::to_string(sub_rows.size() - 1));
            } else {
                // Logger::log(" - Too small, skipped");
            }
        }
        current_x = interval.second;
    }
    
    // Create sub-row after the last blockage (if any space left)
    if (current_x < row_end_x) {
        double sub_start = current_x;
        double sub_end = row_end_x;
        
        // Logger::log("  Trying to create sub-row after last blockage: [" 
                //  + std::to_string(sub_start) + ", " + std::to_string(sub_end) + "]");
        
        // Ensure sub-row boundaries are site-aligned
        // sub_start: ceil to next site boundary
        double relative_start = sub_start - row_start_x;
        int start_site_offset = static_cast<int>(std::ceil(relative_start / site_width));
        sub_start = row_start_x + start_site_offset * site_width;
        
        // sub_end: should already be site-aligned (end_x = start_x + site_width * num_sites)
        // but apply floor for consistency
        double relative_end = sub_end - row_start_x;
        int end_site_offset = static_cast<int>(std::floor(relative_end / site_width));
        sub_end = row_start_x + end_site_offset * site_width;
        
        // Logger::log("    After site alignment: [" 
                //  + std::to_string(sub_start) + ", " + std::to_string(sub_end) + "]");
        
        if (sub_end > sub_start) {
            sub_rows.emplace_back(row_index, sub_start, sub_end, y, height, site_width);
            // Logger::log(" - Created sub-row " + std::to_string(sub_rows.size() - 1));
        } else {
            // Logger::log(" - Too small, skipped");
        }
    }
    
    // Logger::log("  Total sub-rows created: " + std::to_string(sub_rows.size()));
}

// PlacementData methods
void PlacementData::initialize_sub_row_pointers() {
    all_sub_rows.clear();
    for (auto& row : rows) {
        for (auto& sub_row : row.sub_rows) {
            all_sub_rows.push_back(&sub_row);
        }
    }
    
    // Logger::log("Initialized " + std::to_string(all_sub_rows.size()) + " sub-row pointers:");
    for (int i = 0; i < static_cast<int>(all_sub_rows.size()); ++i) {
        // Logger::log("  Sub-row " + std::to_string(i) + ": y=" + std::to_string(all_sub_rows[i]->y) 
                //  + ", x_range=[" + std::to_string(all_sub_rows[i]->start_x) 
                //  + ", " + std::to_string(all_sub_rows[i]->end_x) 
                //  + "], sites=" + std::to_string(all_sub_rows[i]->available_sites));
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