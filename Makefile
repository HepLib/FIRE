default:
	$(MAKE) -C ./poly
	$(MAKE) -C ./prime
	$(MAKE) -C ./MPQ

all:
	$(MAKE) -C ./poly
	$(MAKE) -C ./prime
	$(MAKE) -C ./mpi
	$(MAKE) -C ./MPQ

poly:
	$(MAKE) -C ./poly

prime:
	$(MAKE) -C ./prime

MPQ:
	$(MAKE) -C ./MPQ

.PHONY: mpi poly prime MPQ
mpi:
	$(MAKE) -C ./mpi

test:
	make -f tests.makefile

.PHONY: cleanall
cleanall: clean cleandepall

.PHONY: clean
clean:
	rm -f mpi/*.o
	rm -f poly/*.o
	rm -f prime/*.o
	rm -f MPQ/*.o

dep: depend

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
ASIE = AS_INTEGRATED_ASSEMBLER=1
else
ASIE =
endif

depend:
	cd extra/snappy-1.1.7/ && mkdir -p build && cd build && cmake ../
	$(MAKE) -C ./extra/snappy-1.1.7/build
	$(MAKE) -C ./extra/snappy-1.1.7/build DESTDIR=../../../ install
	#cp extra/kyotocabinet-1.2.77/kccachedb.h.threadsafe.in extra/kyotocabinet-1.2.77/kccachedb.h
	cp extra/kyotocabinet-1.2.77/kccachedb.h.in extra/kyotocabinet-1.2.77/kccachedb.h
	cd extra/kyotocabinet-1.2.77/ && ./configure --prefix=`pwd`/../../usr/
	$(ASIE) $(MAKE) -C ./extra/kyotocabinet-1.2.77
	$(MAKE) -C ./extra/kyotocabinet-1.2.77 install
	$(MAKE) -C ./extra/lz4-1.9.2
	$(MAKE) -C ./extra/lz4-1.9.2 PREFIX=`pwd`/usr/ uninstall
	$(MAKE) -C ./extra/lz4-1.9.2 PREFIX=`pwd`/usr/ install
	cd extra/gperftools-2.7/ && ./configure --enable-minimal --prefix=`pwd`/../../usr/
	$(MAKE) -C ./extra/gperftools-2.7
	$(MAKE) -C ./extra/gperftools-2.7 install
	cp extra/FSBAllocator.hh usr/include/
	cd extra/gmp-6.2.1 && ./configure --enable-cxx --prefix=`pwd`/../../usr/
	$(MAKE) -C ./extra/gmp-6.2.1 install
	#$(MAKE) -C ./extra/zstd-1.4.3
	#$(MAKE) -C ./extra/zstd-1.4.3 PREFIX=`pwd`/usr/ uninstall
	#$(MAKE) -C ./extra/zstd-1.4.3 PREFIX=`pwd`/usr/ install

cleandepall: cleandep
	rm -rf usr/include/*
	rm -rf usr/lib/*
	rm -rf usr/share/*
	rm -rf usr/lib64/*
	rm -rf usr/bin/*
	
cleandep:
	$(MAKE) -C ./extra/gperftools-2.7 clean
	if [ -f extra/kyotocabinet-1.2.77/Makefile ]; then $(MAKE) -C ./extra/kyotocabinet-1.2.77 distclean; fi;
	rm -rf extra/snappy-1.1.7/build
	$(MAKE) -C ./extra/lz4-1.9.2 PREFIX=`pwd`/usr/ clean
	$(MAKE) -C ./extra/gmp-6.2.1 clean
