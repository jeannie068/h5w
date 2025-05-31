#ifndef PARSER_HPP
#define PARSER_HPP

#include "../data_structure/data_structure.hpp"
#include <string>
#include <fstream>
#include <iostream>

class Parser {
public:
    // Parse input file and populate PlacementData
    static bool parse_input_file(const std::string& filename, PlacementData& data);
    
    // Write output file with legalization results
    static bool write_output_file(const std::string& filename, const PlacementData& data);

private:
    // Helper methods for parsing
    static bool parse_max_displacement_constraint(std::ifstream& file, PlacementData& data);
    static bool parse_cells(std::ifstream& file, PlacementData& data);
    static bool parse_blockages(std::ifstream& file, PlacementData& data);
    static bool parse_rows(std::ifstream& file, PlacementData& data);
    
    // Helper method to skip empty lines and comments
    static bool get_next_valid_line(std::ifstream& file, std::string& line);
    
    // Helper method to trim whitespace
    static std::string trim(const std::string& str);
    
    // Helper method to split string by whitespace
    static std::vector<std::string> split(const std::string& str);
};

#endif // PARSER_HPP