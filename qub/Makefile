# Compiler flags
MEX = mex
MEXFLAGS = -O -Lqubtree -lqubtree -DMEX_FCN
#-O2

#MEX_EXT = mexglx  #32-bit linux
MEX_EXT = mexa64  #64-bit linux

# Info for QuB source tree, used by qub_milOptimize.cpp
QUBDIR = ../../qub
QUBI = -I${QUBDIR}/include -I${QUBDIR}/_qubtree -I${QUBDIR}/qublib -I${QUBDIR}/qubopt
QUBL = -L${QUBDIR}/QuBSuite

# General dependancies
DEPS = qubmatlab.h qubtree/libqubtree.a

# Primary target: make all MEX functions
TARGETS = qub_loadTree.${MEX_EXT} qub_saveTree.${MEX_EXT} qub_milOptimize.${MEX_EXT}

default: ${TARGETS}

qub_loadTree.${MEX_EXT}: qub_loadTree.cpp treestruct.cpp ${DEPS}
	${MEX} ${MEXFLAGS} qub_loadTree.cpp treestruct.cpp
	
qub_saveTree.${MEX_EXT}: qub_saveTree.cpp treestruct.cpp ${DEPS}
	${MEX} ${MEXFLAGS} qub_saveTree.cpp treestruct.cpp

qub_milOptimize.${MEX_EXT}: qub_milOptimize.cpp treestruct.cpp ${DEPS}
	${MEX} ${MEXFLAGS} qub_milOptimize.cpp treestruct.cpp

clean:
	rm -f ${TARGETS}

