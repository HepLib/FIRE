# Note: need to set USR_DIR and MALLOC

vpath %.cpp ../src
vpath %.h ../src

CC:=g++ -std=c++17 -fopenmp
cc:=gcc
CFLAGS:=-O3 -Wno-unused-result -Wno-deprecated-declarations -I${USR_DIR}/include -DFloatR

LFLAGS:=-Wl,-rpath,${USR_DIR}/lib -L${USR_DIR}/lib/ -lpthread -lm -lc -lz -lgmp -lmpfr -lflint $(MALLOC)
UNAME_M:=$(shell uname -m)
ifeq ($(UNAME_M),arm64)
CFLAGS:=$(CFLAGS) -I/opt/homebrew/include
LFLAGS:=$(LFLAGS) -Wl,-rpath,/opt/homebrew/lib -L/opt/homebrew/lib
endif
LFLAGS:=$(LFLAGS) -lpthread -lm -lc -lz -lgmp -lmpfr -lflint $(MALLOC)

UNAME_S:=$(shell uname -s)
ifeq ($(UNAME_S),Darwin)
LFLAGS:=$(LFLAGS) -lc++
else
LFLAGS:=$(LFLAGS) -lstdc++ -lrt
endif

OBJ = point.o equation.o functions.o gateToFermat.o common.o parser.o primes.o

default: FIRE

FIRE: $(OBJ) main.o
	$(CC) $(OBJ) main.o $(LFLAGS) -o FIRE

.cpp.o:
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f ./*.o

point.o: common.h point.h point.cpp
gateToFermat.o: gateToFermat.h gateToFermat.cpp
equation.o: common.h gateToFermat.h point.h equation.h equation.cpp
main.o: common.h gateToFermat.h point.h equation.h functions.h main.cpp
functions.o: common.h parser.h gateToFermat.h point.h equation.h functions.h functions.cpp
parser.o: common.h point.h equation.h parser.h parser.cpp
