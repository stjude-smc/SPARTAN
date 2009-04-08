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
TARGET = ${TARGET_DIR}/medianfilterx.${MEX_EXT}


default: ${TARGET}

${TARGET}: mfilt.cpp mf_mex.cpp ${DEPS}
	${MEX} ${MEXFLAGS} mfilt.cpp mf_mex.cpp -output ${TARGET}

test: mfilt.cpp mf_mex.cpp test.cpp ${DEPS}
	g++ mfilt.cpp mf_mex.cpp test.cpp -o test

clean:
	rm -f ${TARGET} test




