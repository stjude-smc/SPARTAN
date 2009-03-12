MEX = mex
MEXFLAGS = -g -Lqubtree -lqubtree
#-O2

DEPS = qubmatlab.h qubtree/libqubtree.a

default: qub_loadTree.mexa64 qub_saveTree.mexa64

qub_loadTree.mexa64: qub_loadTree.cpp treestruct.cpp ${DEPS}
	${MEX} ${MEXFLAGS} qub_loadTree.cpp treestruct.cpp
	
qub_saveTree.mexa64: qub_saveTree.cpp treestruct.cpp ${DEPS}
	${MEX} ${MEXFLAGS} qub_saveTree.cpp treestruct.cpp

clean:
	rm -f qub_loadTree.mexa64 qub_saveTree.mexa64

