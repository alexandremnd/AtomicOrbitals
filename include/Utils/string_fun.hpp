#pragma once

#include <string>

/**
 * @brief Converts a string to lowercase
 *
 * @param str String to be lowercased
 * @return std::string Lowercased string
 */
std::string to_lowercase(const std::string &str);

/**
 * @brief Replaces a letter by another one
 *
 */
void replace_letter(std::string &str, char old_char, char new_char);