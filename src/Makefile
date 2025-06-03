CXX      := g++
CXXFLAGS := -std=c++17 -g -O0 -Wall -Wextra -MMD
LIBS     := -lm
EXEC     := ../bin/hw5
SRC_DIRS := .\
            parser\
            data_structure\
            legalizer

SRCS     := $(wildcard $(SRC_DIRS:=/*.cpp))
OBJS     := $(SRCS:.cpp=.o)
DEPS     := $(OBJS:.o=.d)

all: $(EXEC)

$(EXEC): $(OBJS)
	@mkdir -p ../bin
	$(CXX) -o $@ $^ $(LIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(EXEC) $(OBJS) $(DEPS)

debug: CXXFLAGS += -DDEBUG
debug: $(EXEC)

ifeq (test, $(firstword $(MAKECMDGOALS)))
  TESTCASE := $(word 2, $(MAKECMDGOALS))
  $(eval $(TESTCASE):;@:)
endif

test: $(EXEC)
	@echo Testing with $(TESTCASE) 
	./$(EXEC) ../testcase/$(TESTCASE).txt ../output/$(TESTCASE).out 

.PHONY: all clean debug test
-include $(DEPS)