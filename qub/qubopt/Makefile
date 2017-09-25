# Compile qubopt library for qub_loadTree and qub_saveTree mex functions, and
# the miltreeiface executable.

CCP = g++
CCPFLAGS = -fPIC -I../qubtree -DMAXILL_EXPORTS -O3 -DNDEBUG -fno-strict-aliasing
LN = g++
LNFLAGS = -shared -fPIC -lm -lgsl -lgslcblas -g


# To compile with MATLAB 2015a, need a specific gcc version.
UNAME := $(shell uname -s)
ifeq ($(UNAME),Linux)
	CPP = g++-4.7
	LN = g++-4.7
endif


OBJ_DIR = obj
OUT_DIR = ../qubsuite

SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:%.cpp=$(OBJ_DIR)/%.o)


default: $(OUT_DIR)/libqubopt.a

$(OUT_DIR)/libqubopt.a: $(OBJS)
	ar rcs $@ $(OBJS)

$(OUT_DIR)/libqubopt.so: $(OBJS)
	$(LN) $(LNFLAGS) -o $@ $(OBJS)

$(OBJ_DIR)/%.o: %.cpp
	$(CCP) $(CCPFLAGS) -c -o $@ $^ 

clean:
	rm -f $(OBJ_DIR)/* $(OUT_DIR)/libqubopt.a $(OUT_DIR)/libqubopt.so



