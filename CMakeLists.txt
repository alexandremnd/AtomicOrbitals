cmake_minimum_required(VERSION 3.10)
project(AtomicOrbitals)

# Set policy CMP0167 to NEW
cmake_policy(SET CMP0167 NEW)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Include directories
include_directories(include)

# Find Python libraries and includes
find_package(Python3 COMPONENTS Development REQUIRED)
include_directories(${Python3_INCLUDE_DIRS})

# Find Boost libraries
find_package(Boost REQUIRED COMPONENTS system filesystem)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

# Source files
file(GLOB_RECURSE SRC_FILES src/*.cpp)

# Add executable
add_executable(AtomicOrbitals ${SRC_FILES})

# Link libraries
target_link_libraries(AtomicOrbitals ${Boost_LIBRARIES} ${Python3_LIBRARIES})

# Set compile options
target_compile_options(AtomicOrbitals PRIVATE -Wall -Wextra -Wunused-parameter -march=native)

# ===== Testing =====
option(ENABLE_TESTS "Enable building tests" OFF)

if (ENABLE_TESTS)
    enable_testing()
    find_package(GTest REQUIRED)
    include_directories(${GTest_INCLUDE_DIRS})

    file(GLOB_RECURSE TEST_FILES test/*.cpp)
    # Source files
    foreach(TEST_FILE ${TEST_FILES})
        get_filename_component(TEST_NAME ${TEST_FILE} NAME_WE)
        add_executable(${TEST_NAME} ${TEST_FILE})
        target_link_libraries(${TEST_NAME} ${Boost_LIBRARIES} ${Python3_LIBRARIES} GTest::GTest GTest::Main)
        add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
    endforeach()
endif()