# Compiler flags
MEX = mex
MEXFLAGS = -O -DMEX_FCN

MEX_EXT = $(shell mexext)
#MEX_EXT = mexw32  #32-bit Windows
#MEX_EXT = mexglx  #32-bit linux
#MEX_EXT = mexa64  #64-bit linux

# General dependancies
DEPS = 

# Primary target: make all MEX functions
TARGET_DIR = ../../cascade_binary
TARGET = ${TARGET_DIR}/forward_viterbi.${MEX_EXT}


default: ${TARGET}

${TARGET}: forward_viterbi.cpp ${DEPS}
	${MEX} ${MEXFLAGS} forward_viterbi.cpp -output ${TARGET}

clean:
	rm -f ${TARGET}





