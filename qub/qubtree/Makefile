# Compiles qubtree library for linking to qub_loadTree and qubSaveTree mex files
# and the miltreeiface executable

CCP = g++
CCPFLAGS = -fPIC -Wno-deprecated -O2
LN = g++
LNFLAGS = -shared -fPIC -lpthread -g

# To compile with MATLAB 2015a, need a specific gcc version.
UNAME := $(shell uname -s)
ifeq ($(UNAME),Linux)
	CPP = g++-4.7
	LN = g++-4.7
endif

OBJ_DIR = obj
OUT_DIR= ../qubsuite

SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:%.cpp=$(OBJ_DIR)/%.o)


default: $(OUT_DIR)/libqubtree.a

$(OUT_DIR)/libqubtree.a: $(OBJS)
	ar rcs $@ $(OBJS)

$(OUT_DIR)/libqubtree.so: $(OBJS)
	$(LN) $(LNFLAGS) -o $@ $(OBJS)

$(OBJ_DIR)/%.o: %.cpp
	$(CCP) $(CCPFLAGS) -c -o $@ $^

clean:
	rm -f $(OBJ_DIR)/* $(OUT_DIR)/libqubtree.a $(OUT_DIR)/libqubtree.so

