#include <string>
#include <algorithm>
#include <cctype>

#include "Utils//to_lowercase.hpp"

std::string to_lowercase(const std::string &str) {
    std::string lower_str = str;
    std::transform(lower_str.begin(), lower_str.end(), lower_str.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return lower_str;
}