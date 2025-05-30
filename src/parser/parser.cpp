#include "parser.hpp"
#include <sstream>
#include <algorithm>

bool Parser::parse_input_file(const std::string& filename, PlacementData& data) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open input file " << filename << std::endl;
        return false;
    }
    
    // Parse max displacement constraint
    if (!parse_max_displacement_constraint(file, data)) {
        std::cerr << "Error: Failed to parse max displacement constraint" << std::endl;
        return false;
    }
    
    // Parse cells
    if (!parse_cells(file, data)) {
        std::cerr << "Error: Failed to parse cells" << std::endl;
        return false;
    }
    
    // Parse blockages
    if (!parse_blockages(file, data)) {
        std::cerr << "Error: Failed to parse blockages" << std::endl;
        return false;
    }
    
    // Parse rows
    if (!parse_rows(file, data)) {
        std::cerr << "Error: Failed to parse rows" << std::endl;
        return false;
    }
    
    file.close();
    
    // Create sub-rows for each row based on blockages
    for (int i = 0; i < static_cast<int>(data.rows.size()); ++i) {
        data.rows[i].create_sub_rows(data.blockages, i);
    }
    
    // Initialize sub-row pointers
    data.initialize_sub_row_pointers();
    
    return true;
}

bool Parser::write_output_file(const std::string& filename, const PlacementData& data) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open output file " << filename << std::endl;
        return false;
    }
    
    // Write total displacement
    file << "TotalDisplacement " << static_cast<int>(data.calculate_total_displacement()) << std::endl;
    
    // Write max displacement
    file << "MaxDisplacement " << static_cast<int>(data.calculate_max_displacement()) << std::endl;
    
    // Write number of cells
    file << "NumCells " << data.cells.size() << std::endl;
    
    // Write cell positions
    for (const auto& cell : data.cells) {
        file << cell.name << " " << static_cast<int>(cell.current_x) << " " 
             << static_cast<int>(cell.current_y) << std::endl;
    }
    
    file.close();
    return true;
}

bool Parser::parse_max_displacement_constraint(std::ifstream& file, PlacementData& data) {
    std::string line;
    if (!get_next_valid_line(file, line)) {
        return false;
    }
    
    std::vector<std::string> tokens = split(line);
    if (tokens.size() != 2 || tokens[0] != "MaxDisplacementConstraint") {
        return false;
    }
    
    try {
        data.max_displacement_constraint = std::stod(tokens[1]);
    } catch (const std::exception&) {
        return false;
    }
    
    return true;
}

bool Parser::parse_cells(std::ifstream& file, PlacementData& data) {
    std::string line;
    if (!get_next_valid_line(file, line)) {
        return false;
    }
    
    std::vector<std::string> tokens = split(line);
    if (tokens.size() != 2 || tokens[0] != "NumCells") {
        return false;
    }
    
    int num_cells;
    try {
        num_cells = std::stoi(tokens[1]);
    } catch (const std::exception&) {
        return false;
    }
    
    // Parse each cell
    for (int i = 0; i < num_cells; ++i) {
        if (!get_next_valid_line(file, line)) {
            return false;
        }
        
        tokens = split(line);
        if (tokens.size() != 6 || tokens[0] != "Cell") {
            return false;
        }
        
        try {
            std::string name = tokens[1];
            int width = std::stoi(tokens[2]);
            int height = std::stoi(tokens[3]);
            double x = std::stod(tokens[4]);
            double y = std::stod(tokens[5]);
            
            data.cells.emplace_back(name, width, height, x, y);
        } catch (const std::exception&) {
            return false;
        }
    }
    
    return true;
}

bool Parser::parse_blockages(std::ifstream& file, PlacementData& data) {
    std::string line;
    if (!get_next_valid_line(file, line)) {
        return false;
    }
    
    std::vector<std::string> tokens = split(line);
    if (tokens.size() != 2 || tokens[0] != "NumBlockages") {
        return false;
    }
    
    int num_blockages;
    try {
        num_blockages = std::stoi(tokens[1]);
    } catch (const std::exception&) {
        return false;
    }
    
    // Parse each blockage
    for (int i = 0; i < num_blockages; ++i) {
        if (!get_next_valid_line(file, line)) {
            return false;
        }
        
        tokens = split(line);
        if (tokens.size() != 6 || tokens[0] != "Blockage") {
            return false;
        }
        
        try {
            std::string name = tokens[1];
            int width = std::stoi(tokens[2]);
            int height = std::stoi(tokens[3]);
            int x = std::stoi(tokens[4]);
            int y = std::stoi(tokens[5]);
            
            data.blockages.emplace_back(name, width, height, x, y);
        } catch (const std::exception&) {
            return false;
        }
    }
    
    return true;
}

bool Parser::parse_rows(std::ifstream& file, PlacementData& data) {
    std::string line;
    if (!get_next_valid_line(file, line)) {
        return false;
    }
    
    std::vector<std::string> tokens = split(line);
    if (tokens.size() != 2 || tokens[0] != "NumRows") {
        return false;
    }
    
    int num_rows;
    try {
        num_rows = std::stoi(tokens[1]);
    } catch (const std::exception&) {
        return false;
    }
    
    // Parse each row
    for (int i = 0; i < num_rows; ++i) {
        if (!get_next_valid_line(file, line)) {
            return false;
        }
        
        tokens = split(line);
        if (tokens.size() != 7 || tokens[0] != "Row") {
            return false;
        }
        
        try {
            std::string name = tokens[1];
            int site_width = std::stoi(tokens[2]);
            int height = std::stoi(tokens[3]);
            int x = std::stoi(tokens[4]);
            int y = std::stoi(tokens[5]);
            int num_sites = std::stoi(tokens[6]);
            
            data.rows.emplace_back(name, site_width, height, x, y, num_sites);
        } catch (const std::exception&) {
            return false;
        }
    }
    
    return true;
}

bool Parser::get_next_valid_line(std::ifstream& file, std::string& line) {
    while (std::getline(file, line)) {
        line = trim(line);
        // Skip empty lines and comments (lines starting with //)
        if (!line.empty() && line.substr(0, 2) != "//") {
            return true;
        }
    }
    return false;
}

std::string Parser::trim(const std::string& str) {
    size_t start = str.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    
    size_t end = str.find_last_not_of(" \t\r\n");
    return str.substr(start, end - start + 1);
}

std::vector<std::string> Parser::split(const std::string& str) {
    std::vector<std::string> tokens;
    std::istringstream iss(str);
    std::string token;
    
    while (iss >> token) {
        tokens.push_back(token);
    }
    
    return tokens;
}