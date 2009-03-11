MEX = mex
MEXFLAGS = -g -Lqubtree -lqubtree
#-O2


default: qub_loadTree.mexa64 qub_saveTree.mexa64

qub_loadTree.mexa64: qub_loadTree.cpp
	${MEX} ${MEXFLAGS} qub_loadTree.cpp
	
qub_saveTree.mexa64: qub_saveTree.cpp
	${MEX} ${MEXFLAGS} qub_saveTree.cpp

clean:
	rm -f qub_loadTree.mexa64 qub_saveTree.mexa64

