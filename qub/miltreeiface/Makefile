# Compiles the miltreeiface executable interface to QuB's MIL rate optimizer.
# Requires qubtree and qubopt libraries.

CCP = g++
CCPFLAGS = -I../qubopt -I../qubtree -fno-strict-aliasing -Wno-deprecated -fPIC
LN = g++
LNFLAGS = -L../qubsuite -lqubtree -lqubopt -lgsl -lgslcblas -lm -lpthread -g

OUT_DIR = ../qubsuite
LIBRARIES = $(OUT_DIR)/libqubtree.a $(OUT_DIR)/libqubtree.a
TARGET = $(OUT_DIR)/miltreeiface


# To compile with MATLAB 2015a, need a specific gcc version.
UNAME := $(shell uname -s)
ifeq ($(UNAME),Linux)
	CPP = g++-4.7
	LN = g++-4.7
	TARGET = $(OUT_DIR)/miltreeiface_linux
endif
ifeq ($(UNAME),Darwin)
	TARGET = $(OUT_DIR)/miltreeiface_mac
endif


all: ${TARGET}
default: ${TARGET}

${TARGET}: miltreeiface.o ${LIBRARIES}
	${LN} -o $@ miltreeiface.o ${LIBRARIES} ${LNFLAGS} 

%.o: %.cpp
	${CCP} ${CCPFLAGS} -c -o $@ $<

clean:
	rm -f *.o ${TARGET}







