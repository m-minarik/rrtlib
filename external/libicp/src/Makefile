CXX      := -g++ -std=c++17 -O3 -funroll-loops
# CXXFLAGS := -pedantic-errors -Wall -Wextra -Werror -Wno-deprecated-copy -Wno-ignored-attributes
LDFLAGS  := -L/usr/lib -lstdc++ -lm

SRC      := $(wildcard *.cpp)
OBJECTS  := $(filter-out demo.o test.o, $(SRC:%.cpp=%.o))

libICP.a: $(OBJECTS)
	ar ruv libICP.a $(OBJECTS)

all: test library

%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@ $(LDFLAGS)

library: libICP.a

test: test.o $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(LIBS) -o test $^ $(LDFLAGS)

demo: demo.o $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(LIBS) -o demo $^ $(LDFLAGS)

.PHONY: all build clean debug release

debug: CXXFLAGS += -DDEBUG -g
debug: all

release: CXXFLAGS += -O3
release: all
