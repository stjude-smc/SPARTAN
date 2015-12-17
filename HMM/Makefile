# Compiler flags
MEX = mex

# OPTIMIZED VERSION
MEXFLAGS = -O -DMEX_FCN -largeArrayDims

# DEBUG VERSION
#MEXFLAGS = -g -DMEX_FCN -largeArrayDims


#MEX_EXT = $(shell mexext)
#MEX_EXT = mexw32  #32-bit Windows
#MEX_EXT = mexglx  #32-bit linux
MEX_EXT = mexa64  #64-bit linux


# General dependancies
DEPS = 


# Primary target: make all MEX functions
TARGET_DIR = ../binary
TARGET = ${TARGET_DIR}/forward_viterbix.${MEX_EXT}


default: ${TARGET}

${TARGET}: forward_viterbix.cpp ${DEPS}
	${MEX} ${MEXFLAGS} forward_viterbix.cpp -output ${TARGET}

clean:
	rm -f ${TARGET}





