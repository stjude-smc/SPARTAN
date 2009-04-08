CC = gcc
CCFLAGS = 
CCP = g++
CCPFLAGS = -fPIC -I.
#-O2
LN = g++
LNFLAGS = -shared -fPIC -lpthread -g
TARGET_DIR = ../../../cascade_binary/
OBJ = ${TARGET_DIR}/qubtree_obj

OBJS =  ${OBJ}/QUB_QTR_Impl.o \
	${OBJ}/StringsInterned.o \
	${OBJ}/istream_readline.o \
	${OBJ}/QUB_Tree.o \
	${OBJ}/QTR_NodeFunc.o \
	${OBJ}/QUB_Tree_text.o

default: ${TARGET_DIR}/libqubtree.a

libqubtree.a: ${OBJS}
	ar rcs ${TARGET_DIR}/libqubtree.a ${OBJS}

libqubtree.so: ${OBJS}
	${LN} ${LNFLAGS} -o ${TARGET_DIR}/libqubtree.a ${OBJS}

clean:
	rm -f ${OBJ}/*


${OBJ}/QUB_QTR_Impl.o: QUB_QTR_Impl.cpp QUB_QTR_Impl.h StringsInterned.h
	${CCP} ${CCPFLAGS} -c QUB_QTR_Impl.cpp -o ${OBJ}/QUB_QTR_Impl.o
${OBJ}/StringsInterned.o: StringsInterned.cpp StringsInterned.h
	${CCP} ${CCPFLAGS} -c StringsInterned.cpp -o ${OBJ}/StringsInterned.o
${OBJ}/istream_readline.o: istream_readline.cpp istream_readline.h
	${CCP} ${CCPFLAGS} -c istream_readline.cpp -o ${OBJ}/istream_readline.o
${OBJ}/QUB_Tree.o: QUB_Tree.cpp QUB_Tree.h QUB_QTR_Impl.h
	${CCP} ${CCPFLAGS} -c QUB_Tree.cpp -o ${OBJ}/QUB_Tree.o
${OBJ}/QTR_NodeFunc.o: QTR_NodeFunc.cpp QTR_NodeFunc.h QUB_Tree.h QUB_QTR_Impl.h
	${CCP} ${CCPFLAGS} -c QTR_NodeFunc.cpp -o ${OBJ}/QTR_NodeFunc.o
${OBJ}/QUB_Tree_text.o: QUB_Tree_text.cpp QUB_Tree.h QUB_QTR_Impl.h istream_readline.h
	${CCP} ${CCPFLAGS} -c QUB_Tree_text.cpp -o ${OBJ}/QUB_Tree_text.o

