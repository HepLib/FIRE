export Thread=0
export TcMalloc=0

all:
	$(MAKE) -C ./FlintM
	$(MAKE) -C ./FMPQ
	$(MAKE) -C ./MPQ
	$(MAKE) -C ./poly
	$(MAKE) -C ./prime
	$(MAKE) -C ./FlintC
	$(MAKE) -C ./FlintX

cleanall: clean cleanexe cleandepall

clean:
	rm -f poly/*.o
	rm -f prime/*.o
	rm -f MPQ/*.o
	rm -f FMPQ/*.o
	rm -f FlintM/*.o FlintM/fmpz_mpoly_q/*.o
	rm -f FlintC/*.o
	rm -f FlintX/*.o

cleanexe:
	rm -f poly/F*
	rm -f prime/F*
	rm -f MPQ/F*
	rm -f FMPQ/F*
	rm -f FlintM/F*
	rm -f FlintC/F*
	rm -f FlintX/F*

dep: depend

UNAME_S:=$(shell uname -s)
ifeq ($(UNAME_S),Darwin)
ASIE = AS_INTEGRATED_ASSEMBLER=1
else
ASIE =
endif

depend:
	cp FlintM/fmpz_mpoly_q/fmpz_mpoly_q.h usr/include/
	
	cd extra/snappy-1.1.7/ && mkdir -p build && cd build && cmake ../
	$(MAKE) -C ./extra/snappy-1.1.7/build
	$(MAKE) -C ./extra/snappy-1.1.7/build DESTDIR=../../../ install
	#===========
	#cd extra/snappy-1.1.9/ && mkdir -p build && cd build && cmake ../
	#$(MAKE) -C ./extra/snappy-1.1.9/build
	#$(MAKE) -C ./extra/snappy-1.1.9/build DESTDIR=../../../ install
	
#ifeq ($(Thread),1)
#	cp extra/kyotocabinet-1.2.77/kccachedb.h.threadsafe.in extra/kyotocabinet-1.2.77/kccachedb.h
#else
#	cp extra/kyotocabinet-1.2.77/kccachedb.h.in extra/kyotocabinet-1.2.77/kccachedb.h
#endif
	#cd extra/kyotocabinet-1.2.77/ && ./configure --prefix=`pwd`/../../usr/
	#$(ASIE) $(MAKE) -C ./extra/kyotocabinet-1.2.77
	#$(MAKE) -C ./extra/kyotocabinet-1.2.77 install

	cd extra/kyotocabinet-1.2.79/ && ./configure --prefix=`pwd`/../../usr/
	$(ASIE) $(MAKE) -C ./extra/kyotocabinet-1.2.79
	$(MAKE) -C ./extra/kyotocabinet-1.2.79 install
	
	#$(MAKE) -C ./extra/zstd-1.4.3
	#$(MAKE) -C ./extra/zstd-1.4.3 PREFIX=`pwd`/usr/ uninstall
	#$(MAKE) -C ./extra/zstd-1.4.3 PREFIX=`pwd`/usr/ install
	mkdir -p extra/zstd-1.5.2/build/cmake/build ; cd extra/zstd-1.5.2/build/cmake/build && cmake -DCMAKE_INSTALL_PREFIX=`pwd`/../../../../../usr/ .. && $(MAKE) install
	
	$(MAKE) -C ./extra/lz4-1.9.2
	$(MAKE) -C ./extra/lz4-1.9.2 PREFIX=`pwd`/usr/ uninstall
	$(MAKE) -C ./extra/lz4-1.9.2 PREFIX=`pwd`/usr/ install
	
	cd extra/gperftools-2.10/ && ./configure --enable-minimal --prefix=`pwd`/../../usr/
	$(MAKE) -C ./extra/gperftools-2.10
	$(MAKE) -C ./extra/gperftools-2.10 install
		
	# needed by poly version
	cp extra/FSBAllocator.hh usr/include/
	
	cd extra/gmp-6.2.1 && ./configure --enable-cxx --prefix=`pwd`/../../usr/
	$(MAKE) -C ./extra/gmp-6.2.1 install
	cd extra/mpfr-4.1.0 && ./configure --with-gmp=`pwd`/../../usr/ --prefix=`pwd`/../../usr/
	$(MAKE) -C ./extra/mpfr-4.1.0 install
	cd extra/flint-2.9.0 && ./configure --with-gmp=`pwd`/../../usr/ --with-mpfr=`pwd`/../../usr/ --prefix=`pwd`/../../usr/ 
	$(MAKE) -C ./extra/flint-2.9.0 install

cleandepall: cleandep
	rm -rf usr/include/*
	rm -rf usr/lib/*
	rm -rf usr/share/*
	rm -rf usr/lib64
	rm -rf usr/bin/*
	cd usr && ln -s lib lib64
	
cleandep:
	$(MAKE) -C ./extra/gperftools-2.10 clean
	
	if [ -f extra/kyotocabinet-1.2.79/Makefile ]; then $(MAKE) -C ./extra/kyotocabinet-1.2.79 distclean; fi;
	rm -rf extra/snappy-1.1.9/build
	$(MAKE) -C ./extra/lz4-1.9.2 PREFIX=`pwd`/usr/ clean
	rm -rf extra/zstd-1.5.2/build/cmake/build
	
	$(MAKE) -C ./extra/gmp-6.2.1 clean
	$(MAKE) -C ./extra/mpfr-4.1.0 clean
	$(MAKE) -C ./extra/flint-2.9.0 clean
