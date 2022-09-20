vpath %.cpp ../sources
vpath %.h ../sources

CC:=g++
cc:=gcc
CFLAGS:=-O3 -I../usr/include -DWITH_SNAPPY -DWITH_ZLIB -DWITH_ZSTD
ifeq ($(TcMalloc),1)
CFLAGS:=$(CFLAGS) -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free
else ifeq ($(JeMalloc),1)
CFLAGS:=$(CFLAGS) -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free
endif
CFLAGS:=$(CFLAGS) -DFlintM 

LFLAGS:=-Wl,-rpath,$$PWD/../usr/lib -L../usr/lib/ -L../usr/lib64/ -lkyotocabinet -lpthread -lm -lc -lsnappy -lz ../usr/lib/libflint.a -lgmp -lmpfr ../usr/lib/liblz4.a ../usr/lib/libzstd.a
UNAME_S:=$(shell uname -s)
ifeq ($(UNAME_S),Darwin)
LFLAGS:=$(LFLAGS) -lc++
else
LFLAGS:=$(LFLAGS) -lstdc++ -lrt
endif

ifeq ($(TcMalloc),1)
LFLAGS:=$(LFLAGS) -ltcmalloc_minimal
endif

fOBJ = fmpz_mpoly_q/add.o fmpz_mpoly_q/clear.o fmpz_mpoly_q/equal.o fmpz_mpoly_q/is_canonical.o fmpz_mpoly_q/neg.o fmpz_mpoly_q/set_fmpq.o fmpz_mpoly_q/sub_fmpq.o fmpz_mpoly_q/add_fmpq.o fmpz_mpoly_q/div.o fmpz_mpoly_q/init.o fmpz_mpoly_q/mul.o fmpz_mpoly_q/print_pretty.o fmpz_mpoly_q/set_fmpz.o fmpz_mpoly_q/sub_fmpz.o fmpz_mpoly_q/add_fmpz.o fmpz_mpoly_q/div_fmpq.o fmpz_mpoly_q/inlines.o fmpz_mpoly_q/mul_fmpq.o fmpz_mpoly_q/randtest.o fmpz_mpoly_q/set_si.o fmpz_mpoly_q/swap.o fmpz_mpoly_q/canonicalise.o fmpz_mpoly_q/div_fmpz.o fmpz_mpoly_q/inv.o fmpz_mpoly_q/mul_fmpz.o fmpz_mpoly_q/set.o fmpz_mpoly_q/sub.o fmpz_mpoly_q/used_vars.o

fmpz_mpoly_q/%.o : fmpz_mpoly_q/%.c
	$(cc) $(CFLAGS) -o $@ -c $<

OBJ = point.o equation.o functions.o gateToFermat.o common.o parser.o primes.o

default: FIRE6 FLAME6 FTool6

FIRE6: $(OBJ) $(fOBJ) main.o
	$(CC) $(OBJ) $(fOBJ) main.o $(LFLAGS) -o FIRE6

FLAME6: $(OBJ) $(fOBJ) thread.o
	$(CC) $(OBJ) $(fOBJ) thread.o $(LFLAGS) -o FLAME6

FTool6: $(OBJ) $(fOBJ) ftool.o
	$(CC) $(OBJ) $(fOBJ) ftool.o $(LFLAGS) -o FTool6

.cpp.o:
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f ./*.o fmpz_mpoly_q/*.o

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
