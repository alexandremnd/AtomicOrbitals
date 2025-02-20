#include <string>
#include <algorithm>
#include <cctype>

#include "Utils/string_fun.hpp"

std::string to_lowercase(const std::string &str) {
    std::string lower_str = str;
    std::transform(lower_str.begin(), lower_str.end(), lower_str.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return lower_str;
}

void replace_letter(std::string &str, char old_char, char new_char) {
    std::replace(str.begin(), str.end(), old_char, new_char);
}