cmake_minimum_required(VERSION 3.12)
project(AtomicOrbitals)

cmake_policy(SET CMP0167 NEW)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ======== Dependencies ========
include_directories(include)

find_package(Python3 COMPONENTS Development REQUIRED)
include_directories(${Python3_INCLUDE_DIRS})

find_package(Boost REQUIRED COMPONENTS system filesystem)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

# ======== Build ========
file(GLOB_RECURSE SRC_FILES src/*.cpp)
list(FILTER SRC_FILES EXCLUDE REGEX ".*main.cpp$")
add_library(AtomicOrbitalsObjects OBJECT ${SRC_FILES})

add_executable(AtomicOrbitals src/main.cpp $<TARGET_OBJECTS:AtomicOrbitalsObjects>)
target_link_libraries(AtomicOrbitals ${Boost_LIBRARIES} ${Python3_LIBRARIES})
target_compile_options(AtomicOrbitals PRIVATE -w -march=native)

# ======== Tests ========
option(ENABLE_TESTS "Enable building tests" OFF)

if (ENABLE_TESTS)
    enable_testing()
    find_package(GTest REQUIRED)
    include_directories(${GTest_INCLUDE_DIRS})

    file(GLOB_RECURSE TEST_FILES test/*.cpp)
    foreach(TEST_FILE ${TEST_FILES})
        get_filename_component(TEST_NAME ${TEST_FILE} NAME_WE)
        add_executable(${TEST_NAME} ${TEST_FILE} $<TARGET_OBJECTS:AtomicOrbitalsObjects>)
        target_link_libraries(${TEST_NAME} ${Boost_LIBRARIES} ${Python3_LIBRARIES} GTest::GTest GTest::Main)
        add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
    endforeach()
endif()