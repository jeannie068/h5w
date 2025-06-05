#include <iostream>
#include <string>
#include <chrono>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include "parser/parser.hpp"
#include "legalizer/abacus_legalizer.hpp"
#include "Logger.hpp"

// Modified main.cpp section for running legalization
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
    
    // Determine if this is a large dataset
    bool is_large_dataset = data.cells.size() >= 220000;
    if (is_large_dataset) {
        std::cout << "  - Large dataset detected, using optimized settings" << std::endl;
    }
    
    // Run legalization
    std::cout << "\nRunning Abacus legalization..." << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();
    
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
    
    // Check if max displacement constraint is satisfied
    if (max_displacement > data.max_displacement_constraint) {
        std::cout << "  WARNING: Max displacement constraint violated!" << std::endl;
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