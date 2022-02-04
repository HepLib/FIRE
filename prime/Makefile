vpath %.cpp ../sources
vpath %.h ../sources

CC= g++
cc= gcc
CFLAGS = -O3 -Wall -Wno-strict-aliasing -std=c++11 -I../usr/include -DPRIME -DWITH_SNAPPY -DWITH_ZLIB
#CFLAGS = -O3 -Wall -Wno-strict-aliasing -std=c++11 -I../usr/include -DPRIME -DWITH_SNAPPY -DWITH_ZLIB -DFSBALLOCATOR_USE_THREAD_SAFE_LOCKING_PTHREAD -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
LFLAGS = -L../usr/lib/ -L../usr/lib64/ -lkyotocabinet ../usr/lib/liblz4.a -lc++ -lpthread -lm -lc -lsnappy -lz
else
LFLAGS = -L../usr/lib/ -L../usr/lib64/ -Wl,-Bstatic -lkyotocabinet -l:liblz4.a -Wl,-Bdynamic -lstdc++ -lrt -lpthread -lm -lc -lsnappy -lz
endif



OBJ = point.o equation.o functions.o gateToFermat.o common.o parser.o primes.o

default: FIRE6 FLAME6 FTool6

FIRE6: $(OBJ) main.o
	${CC} $(OBJ) main.o ${LFLAGS} -o FIRE6

FLAME6: $(OBJ) thread.o
	${CC} $(OBJ) thread.o ${LFLAGS} -o FLAME6

FTool6: $(OBJ) ftool.o
	${CC} $(OBJ) ftool.o ${LFLAGS} -o FTool6


.cpp.o:
	${CC} ${CFLAGS} ${INCDIR} -c $<

clean:
	rm -f ./*.o test

point.o: point.cpp
gateToFermat.o: gateToFermat.cpp
equation.o: equation.cpp
functions.o: functions.cpp
main.o: main.cpp
point.o: point.h
gateToFermat.o: gateToFermat.h
equation.o: equation.h
equation.o: point.h
equation.o: gateToFermat.h
functions.o: functions.h
functions.o: equation.h
functions.o: point.h
functions.o: gateToFermat.h
main.o: functions.h
main.o: equation.h
main.o: point.h
main.o: gateToFermat.h
point.o: common.h
common.o: common.h
equation.o: common.h
main.o: common.h
functions.o: parser.h
parser.o: common.h
parser.o: point.h
parser.o: equation.h
thread.o: functions.h
thread.o: equation.h
thread.o: point.h
thread.o: gateToFermat.h
ftool.o: functions.h
ftool.o: equation.h
ftool.o: point.h
ftool.o: gateToFermat.h
main.o: handler.h
thread.o: handler.h
