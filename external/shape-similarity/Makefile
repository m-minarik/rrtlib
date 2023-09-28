# Local
CXX      := -g++ -std=c++17 -O3
CXXFLAGS := -Wno-unused-result -Wno-format -Wno-return-type
# CXXFLAGS := -pedantic-errors -Wall -Wextra -Werror -Wno-deprecated-copy -Wno-ignored-attributes

# Metacentrum
# CXX		 := icpc -std=c++17 -O2

INCLUDE  := -Iinclude/
LDFLAGS  := -L/usr/lib -lstdc++ -lm
BUILD	 := ./build
OBJ_DIR  := $(BUILD)/objects
APP_DIR  := $(BUILD)/apps

SRC      := $(wildcard src/*.cpp)
OBJECTS  := $(filter-out $(OBJ_DIR)/src/Main.o $(OBJ_DIR)/src/identify_object.o,\
						 $(SRC:%.cpp=$(OBJ_DIR)/%.o))

all: build main identify-object

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@ $(LDFLAGS)

main: $(OBJ_DIR)/src/Main.o $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -o main $^ $(LDFLAGS)

identify-object: $(OBJ_DIR)/src/identify_object.o $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -o identify-object $^ $(LDFLAGS)


.PHONY: all build clean debug release

build:
	@mkdir -p $(OBJ_DIR)

debug: CXXFLAGS += -DDEBUG -g
debug: all

release: CXXFLAGS += -O2
release: all

clean:
	@rm -rvf $(OBJ_DIR)/*
