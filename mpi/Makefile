vpath %.cpp ../sources

TARGETS=FIRE_MPI
MCXX=mpicxx
CXXFLAGS=-O3 -std=c++11 -I../usr/include -DPRIME -Wall

$(TARGETS): mpi_wrapper.o
	${MCXX} -o $@ $<

%.o: %.cpp
	${MCXX} ${CXXFLAGS} ${INCDIR} -c $<

.PHONY: clean
clean:
	rm -f ./*.o $(TARGETS)

mpi_wrapper.o: Makefile
