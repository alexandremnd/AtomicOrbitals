cmake_minimum_required(VERSION 3.10)
project(AtomicOrbitals)

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