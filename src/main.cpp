#include <iostream>
#include <string>
#include <chrono>
#include <algorithm>
#include <cmath>
#include "parser/parser.hpp"
#include "legalizer/abacus_legalizer.hpp"
#include "Logger.hpp"

int main(int argc, char* argv[]) {
    // Check command line arguments
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <input_file> <output_file>" << std::endl;
        std::cerr << "Example: " << argv[0] << " ../testcase/public1.txt ../output/public1.out" << std::endl;
        return 1;
    }
    
    std::string input_file = argv[1];
    std::string output_file = argv[2];
    
    // Parse input file
    PlacementData data;
    std::cout << "Reading input file: " << input_file << std::endl;
    
    if (!Parser::parse_input_file(input_file, data)) {
        std::cerr << "Error: Failed to parse input file" << std::endl;
        return 1;
    }
    
    std::cout << "Input file parsed successfully:" << std::endl;
    std::cout << "  - Number of cells: " << data.cells.size() << std::endl;
    std::cout << "  - Number of blockages: " << data.blockages.size() << std::endl;
    std::cout << "  - Number of rows: " << data.rows.size() << std::endl;
    std::cout << "  - Max displacement constraint: " << data.max_displacement_constraint << std::endl;
    
    // Count total sub-rows
    int total_sub_rows = 0;
    for (const auto& row : data.rows) {
        total_sub_rows += row.sub_rows.size();
    }
    std::cout << "  - Total sub-rows (after blockage handling): " << total_sub_rows << std::endl;
    
    // Run legalization
    std::cout << "\nRunning Abacus legalization..." << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();

    Logger::init();
    
    AbacusLegalizer legalizer;
    bool success = legalizer.legalize(data);
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    if (!success) {
        std::cerr << "Error: Legalization failed" << std::endl;
        return 1;
    }
    
    std::cout << "Legalization completed successfully!" << std::endl;
    std::cout << "  - Runtime: " << duration.count() / 1000.0 << " seconds" << std::endl;
    
    // Calculate and display results
    double total_displacement = data.calculate_total_displacement();
    double max_displacement = data.calculate_max_displacement();
    
    std::cout << "\nResults:" << std::endl;
    std::cout << "  - Total displacement: " << static_cast<int>(total_displacement) << std::endl;
    std::cout << "  - Max displacement: " << static_cast<int>(max_displacement) << std::endl;
    
    // Verify constraints
    bool constraints_satisfied = true;
    
    // Check max displacement constraint
    for (const auto& cell : data.cells) {
        if (cell.get_displacement() > data.max_displacement_constraint) {
            std::cerr << "Warning: Cell " << cell.name << " violates max displacement constraint: "
                     << cell.get_displacement() << " > " << data.max_displacement_constraint << std::endl;
            constraints_satisfied = false;
        }
    }
    
    // Check non-overlapping constraint
    for (const auto& row : data.rows) {
        for (const auto& sub_row : row.sub_rows) {
            std::vector<int> cell_indices = sub_row.cells;
            
            // Sort cells by current x position
            std::sort(cell_indices.begin(), cell_indices.end(),
                      [&data](int a, int b) {
                          return data.cells[a].current_x < data.cells[b].current_x;
                      });
            
            for (int i = 1; i < static_cast<int>(cell_indices.size()); ++i) {
                int prev_idx = cell_indices[i-1];
                int curr_idx = cell_indices[i];
                
                double prev_end = data.cells[prev_idx].current_x + data.cells[prev_idx].width;
                if (prev_end > data.cells[curr_idx].current_x + 1e-9) {
                    std::cerr << "Warning: Cells " << data.cells[prev_idx].name 
                             << " and " << data.cells[curr_idx].name << " overlap!" << std::endl;
                    constraints_satisfied = false;
                }
            }
        }
    }
    
    // Check site alignment constraint
    for (const auto& cell : data.cells) {
        // Find which sub-row contains this cell
        for (const auto& row : data.rows) {
            for (const auto& sub_row : row.sub_rows) {
                if (std::find(sub_row.cells.begin(), sub_row.cells.end(), 
                             &cell - &data.cells[0]) != sub_row.cells.end()) {
                    // Check if cell is site-aligned
                    double relative_x = cell.current_x - sub_row.start_x;
                    int site_offset = static_cast<int>(std::round(relative_x / sub_row.site_width));
                    double aligned_x = sub_row.start_x + site_offset * sub_row.site_width;
                    
                    if (std::abs(cell.current_x - aligned_x) > 1e-9) {
                        std::cerr << "Warning: Cell " << cell.name 
                                 << " is not site-aligned!" << std::endl;
                        constraints_satisfied = false;
                    }
                    break;
                }
            }
        }
    }
    
    if (constraints_satisfied) {
        std::cout << "All constraints satisfied!" << std::endl;
    } else {
        std::cerr << "Warning: Some constraints are violated!" << std::endl;
    }
    
    // Write output file
    std::cout << "\nWriting output file: " << output_file << std::endl;
    if (!Parser::write_output_file(output_file, data)) {
        std::cerr << "Error: Failed to write output file" << std::endl;
        return 1;
    }
    
    std::cout << "Output file written successfully!" << std::endl;
    
    return 0;
}