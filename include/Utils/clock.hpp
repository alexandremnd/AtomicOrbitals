#pragma once

#include <chrono>
#include <iostream>

class Clock {
public:
    Clock() : m_start(std::chrono::high_resolution_clock::now()) {}

    void reset() {
        m_start = std::chrono::high_resolution_clock::now();
    }

    void time(const std::string& message) {
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - m_start);
        std::cout << message << " " << duration.count() << "us" << std::endl;
    }
private:
    std::chrono::time_point<std::chrono::high_resolution_clock> m_start;
};