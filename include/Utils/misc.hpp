#pragma once

#include <memory>
#include <string>
#include <fstream>

std::unique_ptr<std::ofstream> create_output_file(std::string output,
                                                  int precision = 10);

void execute_python_script(std::string name);