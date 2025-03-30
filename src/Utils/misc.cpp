#include "Utils/misc.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

namespace fs = std::filesystem;

std::unique_ptr<std::ofstream> create_output_file(std::string output,
                                                  int precision) {
    fs::path path = fs::current_path() / "data" / output;
    fs::create_directories(path.parent_path());

    auto output_file = std::make_unique<std::ofstream>(path);
    if (!output_file->is_open()) {
        std::cerr << "Output file creation: Could not open the output file "
                  << path << ".\n";
        exit(1);
    }
    *output_file << std::setprecision(precision);

    return output_file;
}

void execute_python_script(std::string name) {
    fs::path path = fs::current_path() / "scripts" / name;
    std::string command = "python3 " + path.string();

    int result = system(command.c_str());

    if (result != 0) {
        std::cerr << "Error executing Python script: " << name << "\n";
        exit(1);
    }
}