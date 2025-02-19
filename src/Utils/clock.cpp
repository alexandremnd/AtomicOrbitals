#include "Utils/clock.hpp"
#include <chrono>
#include <iostream>

void Clock::start() { m_start = std::chrono::high_resolution_clock::now(); }

void Clock::reset() { m_start = std::chrono::high_resolution_clock::now(); }

void Clock::time_us(const std::string &message) {
    auto end = std::chrono::high_resolution_clock::now();
    auto duration =
        std::chrono::duration_cast<std::chrono::microseconds>(end - m_start);
    std::cout << message << " " << duration.count() << "us" << std::endl;
}

void Clock::time_ms(const std::string &message) {
    auto end = std::chrono::high_resolution_clock::now();
    auto duration =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - m_start);
    std::cout << message << " " << duration.count() << "ms" << std::endl;
}

void Clock::time_s(const std::string &message) {
    auto end = std::chrono::high_resolution_clock::now();
    auto duration =
        std::chrono::duration_cast<std::chrono::seconds>(end - m_start);
    std::cout << message << " " << duration.count() << "s" << std::endl;
}