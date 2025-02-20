cmake_minimum_required(VERSION 3.12)
project(AtomicOrbitals)

# cmake_policy(SET CMP0167 NEW)
option(ENABLE_TESTS "Enable building tests" OFF)
option(ENABLE_DOCS "Enable building documentation" OFF)

# ======== Configuration ========

# Tests require a Debug build type to keep assertions enabled
if (ENABLE_TESTS)
    set(CMAKE_BUILD_TYPE Debug)
endif()


set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -fopenmp")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")

# ======== Dependencies ========
include_directories(include)

find_package(Python3 COMPONENTS Development REQUIRED)
include_directories(${Python3_INCLUDE_DIRS})

find_package(Boost REQUIRED COMPONENTS system filesystem)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

find_package(OpenMP)
if (OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
endif()

# ======== Build ========
file(GLOB_RECURSE SRC_FILES "${PROJECT_SOURCE_DIR}/src/*.cpp")
list(FILTER SRC_FILES EXCLUDE REGEX ".*main.cpp$")
add_library(AtomicOrbitalsObjects OBJECT ${SRC_FILES})

add_executable(AtomicOrbitals "${PROJECT_SOURCE_DIR}/src/main.cpp" $<TARGET_OBJECTS:AtomicOrbitalsObjects>)
target_link_libraries(AtomicOrbitals ${Boost_LIBRARIES} ${Python3_LIBRARIES})
target_compile_options(AtomicOrbitals PRIVATE -w)

# ======== Tests ========
if (ENABLE_TESTS)
    enable_testing()
    find_package(GTest REQUIRED)
    include_directories(${GTest_INCLUDE_DIRS})

    file(GLOB_RECURSE TEST_FILES "${PROJECT_SOURCE_DIR}/test/*.cpp")
    foreach(TEST_FILE ${TEST_FILES})
        get_filename_component(TEST_NAME ${TEST_FILE} NAME_WE)
        add_executable(${TEST_NAME} ${TEST_FILE} $<TARGET_OBJECTS:AtomicOrbitalsObjects>)
        target_link_libraries(${TEST_NAME} ${Boost_LIBRARIES} ${Python3_LIBRARIES} GTest::GTest GTest::Main)
        set_target_properties(${TEST_NAME} PROPERTIES COMPILE_FLAGS "-DDEBUG")
        add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
    endforeach()
endif()

# ======== Documentation ========
if (ENABLE_DOCS)
    find_package(Doxygen REQUIRED)

    set(DOXYGEN_GENERATE_HTML YES)
    set(DOXYGEN_CONFIG_FILE ${PROJECT_SOURCE_DIR}/docs/Doxyfile)
    set(DOXYGEN_OUTPUT_CONFIG_FILE ${PROJECT_BINARY_DIR}/Doxyfile)
    set(DOXYGEN_INPUT_DIRS "${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/src")

    configure_file(${DOXYGEN_CONFIG_FILE} ${DOXYGEN_OUTPUT_CONFIG_FILE} @ONLY)

    add_custom_target(docs ALL
        COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen" VERBATIM
    )

endif()

message(STATUS "=========== Configuration ===========")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Tests: ${ENABLE_TESTS}")
message(STATUS "Documentation: ${ENABLE_DOCS}")
message(STATUS "Compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "C++ standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "OpenMP available: ${OpenMP_CXX_FOUND}")
message(STATUS "Source directory: ${PROJECT_SOURCE_DIR}")
message(STATUS "Binary directory: ${PROJECT_BINARY_DIR}")
message(STATUS "=========== End Configuration ===========")