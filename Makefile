# Local
CXX      := -g++ -std=c++17 -O2
# CXXFLAGS := -pedantic-errors -Wall -Wextra -Werror -Wno-deprecated-copy -Wno-ignored-attributes
LDFLAGS  := -L/usr/lib -lstdc++ -lm -Lexternal/rapid -lRAPID -L./external/libicp/src -lICP -lboost_system -lboost_serialization -lboost_filesystem -lompl
INCLUDE  := -Iinclude/ -Iexternal/rapid -Iexternal/nanoflann/include/ -Iexternal/json/single_include/nlohmann/ -Iexternal/libicp/src/

# EIGEN setup
INCLUDE  += -I/usr/local/include/

# OMPL setup
LDFLAGS  += -L$(HOME)/opt/lib 
INCLUDE  += -I$(HOME)/opt/include

BUILD	 := ./build
OBJ_DIR  := $(BUILD)/objects

SRC      := $(wildcard src/*.cpp)
OBJECTS  := $(SRC:%.cpp=$(OBJ_DIR)/%.o)

all: build main

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@ $(LDFLAGS)
	
.PHONY: all build clean debug release main

build:
	@mkdir -p $(OBJ_DIR)

debug: CXXFLAGS += -DDEBUG -g
debug: all

release: CXXFLAGS += -O3
release: all

clean:
	@rm -rvf $(OBJ_DIR)/*

main: $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(LIBS) -o main $^ $(LDFLAGS)