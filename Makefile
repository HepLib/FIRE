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

#=====================================================

dep: depend

UNAME_S:=$(shell uname -s)
ifeq ($(UNAME_S),Darwin)
ASIE = AS_INTEGRATED_ASSEMBLER=1
else
ASIE =
endif

USR_DIR:=$(shell pwd)/usr
depend:
	cd extra; \
tar zxf gmp-6.2.1.tar.gz; \
tar zxf mpfr-4.1.1.tar.gz; \
tar zxf flint-2.9.0.tar.gz; \
tar zxf snappy-1.1.7.tar.gz; \
tar zxf zstd-1.5.2.tar.gz; \
tar zxf lz4-1.9.4.tar.gz; \
tar zxf kyotocabinet-1.2.79.tar.gz

	cp FlintM/fmpz_mpoly_q/fmpz_mpoly_q.h usr/include/
	
	cd extra/gmp-6.2.1;\
./configure --enable-cxx --prefix=$(USR_DIR); \
$(MAKE) install
	cd extra/mpfr-4.1.1; \
./configure --with-gmp=$(USR_DIR) --prefix=$(USR_DIR); \
$(MAKE) install
	cd extra/flint-2.9.0; \
./configure --with-gmp=$(USR_DIR) --with-mpfr=$(USR_DIR) --prefix=$(USR_DIR); \
$(MAKE) install
	
	cd extra/snappy-1.1.7; \
mkdir -p build; \
cd build; \
cmake -DCMAKE_INSTALL_PREFIX=$(USR_DIR) .. ; \
$(MAKE) install
	
	cd extra/kyotocabinet-1.2.79; \
./configure --prefix=$(USR_DIR); \
$(ASIE) $(MAKE); \
$(MAKE) install
	
	mkdir -p extra/zstd-1.5.2/build/cmake/build; \
cd extra/zstd-1.5.2/build/cmake/build; \
cmake -DCMAKE_INSTALL_PREFIX=$(USR_DIR) ..; \
$(MAKE) install
	
	$(MAKE) -C ./extra/lz4-1.9.4 PREFIX=$(USR_DIR) install
			
	# needed by poly version
	cp extra/FSBAllocator.hh usr/include/
	
cleandepall: cleandep
	rm -rf usr/include/*
	rm -rf usr/lib/*
	rm -rf usr/share/*
	rm -rf usr/lib64
	rm -rf usr/bin/*
	cd usr && ln -s lib lib64
	
cleandep:
	rm -rf extra/gmp-6.2.1
	rm -rf extra/mpfr-4.1.1
	rm -rf extra/flint-2.9.0
	rm -rf extra/snappy-1.1.7
	rm -rf extra/zstd-1.5.2
	rm -rf extra/lz4-1.9.4
	rm -rf extra/kyotocabinet-1.2.79
