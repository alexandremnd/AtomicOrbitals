#pragma once

#include <chrono>

class Clock {
  public:
    Clock() : m_start(std::chrono::high_resolution_clock::now()) {}

    void start();
    void reset();
    void time_us(const std::string &message);
    void time_ms(const std::string &message);
    void time_s(const std::string &message);

  private:
    std::chrono::time_point<std::chrono::high_resolution_clock> m_start;
};