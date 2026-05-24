CXX = mpicxx
CXXFLAGS = -Wall -Wextra -O2 -std=c++17

CORE_TARGET = solver-core
EXT_TARGET = solver-extension
OBJDIR = .obj

CORE_SRCS = $(wildcard core/*.cpp)
EXT_SRCS = $(wildcard extension/*.cpp)

CORE_OBJS = $(patsubst core/%.cpp,$(OBJDIR)/core/%.o,$(CORE_SRCS))
EXT_OBJS = $(patsubst extension/%.cpp,$(OBJDIR)/extension/%.o,$(EXT_SRCS))

.PHONY: all build build-core build-ext format clean

all: build

build: build-core build-ext

build-core: $(CORE_TARGET)

build-ext: $(EXT_TARGET)

$(CORE_TARGET): $(CORE_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(EXT_TARGET): $(EXT_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(OBJDIR)/core/%.o: core/%.cpp
	@mkdir -p $(OBJDIR)/core
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(OBJDIR)/extension/%.o: extension/%.cpp
	@mkdir -p $(OBJDIR)/extension
	$(CXX) $(CXXFLAGS) -c -o $@ $<

format:
	clang-format -i $(CORE_SRCS) $(EXT_SRCS)

clean:
	rm -f $(CORE_TARGET) $(EXT_TARGET)
	rm -rf $(OBJDIR)
