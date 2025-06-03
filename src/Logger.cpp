// Logger.cpp
#include "Logger.hpp"

// Initialize static members
std::ofstream Logger::logFile;
bool Logger::initialized = false;
int Logger::indent = 0;