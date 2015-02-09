# Compiles the qub_loadTree and qub_saveTree MEX files for loading and saving
# QuB-format files -- specifically qmf model files.
# Requires mex be in the path and will not work with Windows.
# Use qubtree_compile_mex.m instead.

# Uncomment the line that corresponds to your system architecture.
#MEX_EXT = $(shell mexext)
#MEX_EXT = mexw32  #32-bit Windows
MEX_EXT = mexw64  #64-bit Windows
#MEX_EXT = mexglx  #32-bit linux
#MEX_EXT = mexa64  #64-bit linux
#MEX_EXT = mexmaci64  #64-bit Mac


# Primary target: make all MEX functions
TARGET_DIR = ../binary
TARGETS = ${TARGET_DIR}/qub_loadTree.${MEX_EXT} ${TARGET_DIR}/qub_saveTree.${MEX_EXT}


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

