#pragma once

#include <chrono>

/**
 * @brief A simple clock class to time code execution.
 *
 */
class Clock {
  public:
    /**
     * @brief Instiates and starts the clock.
     */
    Clock() : m_start(std::chrono::high_resolution_clock::now()) {}

    /**
     * @brief Starts the clock.
     */
    void start();

    /**
     * @brief Starts again the clock
     */
    void reset();

    /**
     * @brief Shows elapsed time in microseconds
     *
     * @param message Message to show before the time
     */
    void time_us(const std::string &message);

    /**
     * @brief Shows elapsed time in milliseconds
     *
     * @param message Message to show before the time
     */
    void time_ms(const std::string &message);

    /**
     * @brief Shows elapsed time in seconds
     *
     * @param message Message to show before the time
     */
    void time_s(const std::string &message);

  private:
    std::chrono::time_point<std::chrono::high_resolution_clock> m_start;
};