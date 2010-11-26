# THIS VERSION does not depend on QuB directly -- 
# MIL is run using a .exe interface program instead that is directly linked.

# Primary target: make all MEX functions
TARGET_DIR = ../../cascade_binary
TARGETS = ${TARGET_DIR}/qub_loadTree.${MEX_EXT} ${TARGET_DIR}/qub_saveTree.${MEX_EXT}


#MEX_EXT = $(shell mexext)
#MEX_EXT = .mexw32  #32-bit Windows
#MEX_EXT = .mexglx  #32-bit linux
MEX_EXT = mexa64  #64-bit linux

# Compiler flags
MEX = mex
MEXFLAGS = -DMEX_FCN -lqubtree_static -Lqubtree

# General dependancies
DEPS = qubmatlab.h qubtree/libqubtree_static.a

default: ${TARGETS}

${TARGET_DIR}/qub_loadTree.${MEX_EXT}: qub_loadTree.cpp treestruct.cpp ${DEPS}
	${MEX} ${MEXFLAGS} qub_loadTree.cpp treestruct.cpp -o $@
	
${TARGET_DIR}/qub_saveTree.${MEX_EXT}: qub_saveTree.cpp treestruct.cpp ${DEPS}
	${MEX} ${MEXFLAGS} qub_saveTree.cpp treestruct.cpp -o $@

clean:
	rm -f ${TARGETS}

