# THIS VERSION does not depend on QuB directly -- 
# MIL is run using a .exe interface program instead that is directly linked.

# Compiler flags
MEX = mex
MEXFLAGS = -O -L../../cascade_binary -Lqubtree -lqubtree -DMEX_FCN

MEX_EXT = $(shell mexext)
#MEX_EXT = mexw32  #32-bit Windows
#MEX_EXT = mexglx  #32-bit linux
#MEX_EXT = mexa64  #64-bit linux

# Primary target: make all MEX functions
TARGET_DIR = ../../cascade_binary
TARGETS = ${TARGET_DIR}/qub_loadTree.${MEX_EXT} ${TARGET_DIR}/qub_saveTree.${MEX_EXT} ${TARGET_DIR}/qub_milOptimize.${MEX_EXT}

# General dependancies
DEPS = qubmatlab.h ${TARGET_DIR}/libqubtree.a

default: ${TARGETS}

${TARGET_DIR}/qub_loadTree.${MEX_EXT}: qub_loadTree.cpp treestruct.cpp ${DEPS}
	${MEX} ${MEXFLAGS} qub_loadTree.cpp treestruct.cpp
	
${TARGET_DIR}/qub_saveTree.${MEX_EXT}: qub_saveTree.cpp treestruct.cpp ${DEPS}
	${MEX} ${MEXFLAGS} qub_saveTree.cpp treestruct.cpp

${TARGET_DIR}/qub_milOptimize.${MEX_EXT}: qub_milOptimize.cpp treestruct.cpp ${DEPS}
	${MEX} ${MEXFLAGS} qub_milOptimize.cpp treestruct.cpp

clean:
	rm -f ${TARGETS}

