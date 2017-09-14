CCP = g++-4.7
CCPFLAGS = -I../qubopt -I../qubtree -fno-strict-aliasing -Wno-deprecated -fPIC
LN = g++-4.7
LNFLAGS = -L../qubsuite -static -lqubtree -lqubopt -lgsl -lgslcblas -lm -lpthread -pthread -g


OUT_DIR = ../qubsuite
LIBRARIES = $(OUT_DIR)/libqubtree.a $(OUT_DIR)/libqubtree.a
TARGET = $(OUT_DIR)/miltreeiface


all: ${TARGET}
default: ${TARGET}

${TARGET}: miltreeiface.o ${LIBRARIES}
	${LN} -o $@ miltreeiface.o ${LIBRARIES} ${LNFLAGS} 

%.o: %.cpp
	${CCP} ${CCPFLAGS} -c -o $@ $<

clean:
	rm -f *.o ${TARGET}







