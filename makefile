CXX := g++
PYTHON_LIBS := $(shell python3-config --ldflags)
BOOST_LIBS := -lboost_system -lboost_filesystem
PYTHON_INCLUDE := $(shell python3-config --includes)
LIBRARIES := $(BOOST_LIBS) $(PYTHON_LIBS) # -fopenmp -llapack -llapacke -lopenblas
CXXFLAGS := -I./include/ $(PYTHON_INCLUDE) -Wall -Wextra -Wunused-parameter -std=c++2b $(LIBRARIES) -march=native # may add -o3 -ffast-math to benchmark
SRC_DIR := ./src
INC_DIR := ./include
BUILD_DIR := ./build
BIN_DIR := ./
TARGET := ${BIN_DIR}/atomOrbital

SRC_FILES := $(shell find $(SRC_DIR) -name '*.cpp')
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(SRC_FILES))

all: $(TARGET)

$(TARGET): $(OBJ_FILES) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

clean:
	rm -rf $(OBJ_DIR)/*.o $(TARGET)

.PHONY: all clean ab