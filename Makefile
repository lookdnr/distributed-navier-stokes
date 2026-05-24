# Compiler and flags
CXX = mpicxx
CXXFLAGS = -Wall -Wextra -O2 -std=c++17

# Targets
CORE_TARGET = solver-core
EXT_TARGET = solver-extension

# Sources
CORE_SRCS = $(wildcard core/*.cpp)
EXT_SRCS = $(wildcard extension/*.cpp)

# Object files
CORE_OBJS = $(CORE_SRCS:.cpp=.o)
EXT_OBJS = $(EXT_SRCS:.cpp=.o)

# Phony targets
.PHONY: all build build-core build-ext format clean

all: build

## Build both executables
build: build-core build-ext

## Build only core
build-core: $(CORE_TARGET)

## Build only extension
build-ext: $(EXT_TARGET)

$(CORE_TARGET): $(CORE_OBJS)	
	$(CXX) $(CXXFLAGS) -o $@ $

$(EXT_TARGET): $(EXT_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $

## compile command
COMPILE = $(CXX) $(CXXFLAGS) -c -o $@ $(patsubst %.o,%.cpp,$@)

%.o: %.cpp
	$(COMPILE)

## Format all source files
format:
	clang-format -i $(CORE_SRCS) $(EXT_SRCS) $(wildcard core/*.hpp extension/*.hpp)

## Clean build artifacts
clean:
	rm -f $(CORE_OBJS) $(EXT_OBJS) $(CORE_TARGET) $(EXT_TARGET)
