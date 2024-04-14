ifeq ($(strip $(MALLOC)),)
MALLOC:=-ljemalloc
endif

ifeq ($(strip $(USR_DIR)),)
USR_DIR:=$(PWD)/usr
endif

ifeq ($(strip $(USR_DIR)),)
USR_DIR:=$(shell pwd)/usr
endif

export MALLOC USR_DIR

defalut: oM oQ oF oP

oM:
	$(MAKE) -C M
oQ:
	$(MAKE) -C Q
oF:
	$(MAKE) -C F
oP:
	$(MAKE) -C P
		
cleanall: clean cleanexe cleandepall

clean: cleandep
	$(MAKE) clean -C M
	$(MAKE) clean -C Q
	$(MAKE) clean -C F
	$(MAKE) clean -C P

cleanexe:
	rm -f M/FIRE
	rm -f Q/FIRE
	rm -f F/FIRE
	rm -f P/FIRE
	
tgz:
	/bin/rm -rf FIREs FIREs.tar.gz; \
mkdir FIREs; \
cp -r src M Q F P Makefile usr usr-src FIREs/; \
make cleanall -C FIREs; \
/bin/rm -f FIREs/.DS_Store FIREs/*/.DS_Store FIREs/*/*/.DS_Store FIREs/*/*/*/.DS_Store;\
gtar cfz FIREs.tar.gz FIREs; \
/bin/rm -rf FIREs
	
UNAME_S:=$(shell uname -s)
ifeq ($(UNAME_S),Darwin)
ASIE = AS_INTEGRATED_ASSEMBLER=1
else
ASIE =
endif

UNAME_M:=$(shell uname -m)
ifeq ($(UNAME_M),arm64)
#===============================================
dep:
	cd usr-src; \
tar xf flint-3.1.2.tar.gz; \
tar xf jemalloc-5.3.0.tar.bz2; \
tar xf gperftools-2.10.tar.gz; \
tar xf mimalloc-2.1.2.tar.gz;

	cd usr-src/flint-3.1.2; \
./bootstrap.sh; \
./configure --disable-static --with-gmp=/opt/homebrew --with-mpfr=/opt/homebrew --prefix=$(USR_DIR) CFLAGS="-O3"; \
$(MAKE) install

	cd usr-src/jemalloc-5.3.0; \
./configure --prefix=$(USR_DIR); \
$(MAKE) install

	cd usr-src/gperftools-2.10; \
./configure --prefix=$(USR_DIR); \
$(MAKE) install

	cd usr-src/mimalloc-2.1.2; \
mkdir build && cd build; \
cmake -DCMAKE_INSTALL_PREFIX=$(USR_DIR) .. ; \
$(MAKE) install

cleandepall: cleandep
	rm -rf usr/include/*
	rm -rf usr/lib/*
	rm -rf usr/share/*
	rm -rf usr/lib64
	rm -rf usr/bin/*
	cd usr && ln -s lib lib64
	
cleandep:
	rm -rf usr-src/flint-3.1.2
	rm -rf usr-src/jemalloc-5.3.0
	rm -rf usr-src/gperftools-2.10
	rm -rf usr-src/mimalloc-2.1.2
#===============================================
else
#===============================================
dep:
	cd usr-src; \
tar xf gmp-6.2.1.tar.gz; \
tar xf mpfr-4.2.1.tar.gz; \
tar xf flint-3.1.2.tar.gz; \
tar xf jemalloc-5.3.0.tar.bz2; \
tar xf gperftools-2.10.tar.gz; \
tar xf mimalloc-2.1.2.tar.gz;

	cd usr-src/gmp-6.2.1;\
./configure --enable-cxx --prefix=$(USR_DIR); \
$(MAKE) install
	cd usr-src/mpfr-4.2.1; \
./configure --with-gmp=$(USR_DIR) --prefix=$(USR_DIR); \
$(MAKE) install
	cd usr-src/flint-3.1.2; \
./bootstrap.sh; \
./configure --disable-static --enable-avx2 --with-gmp=$(USR_DIR) --with-mpfr=$(USR_DIR) --prefix=$(USR_DIR) CFLAGS="-O3"; \
$(MAKE) install
		
	cd usr-src/jemalloc-5.3.0; \
./configure --prefix=$(USR_DIR); \
$(MAKE) install

	cd usr-src/gperftools-2.10; \
./configure --prefix=$(USR_DIR); \
$(MAKE) install

	cd usr-src/mimalloc-2.1.2; \
mkdir build && cd build; \
cmake -DCMAKE_INSTALL_PREFIX=$(USR_DIR) .. ; \
$(MAKE) install

cleandep:
	rm -rf usr-src/gmp-6.2.1
	rm -rf usr-src/mpfr-4.2.1
	rm -rf usr-src/flint-3.1.2
	rm -rf usr-src/jemalloc-5.3.0
	rm -rf usr-src/gperftools-2.10
	rm -rf usr-src/mimalloc-2.1.2

cleandepall: cleandep
	rm -rf usr/include/*
	rm -rf usr/lib/*
	rm -rf usr/share/*
	rm -rf usr/lib64
	rm -rf usr/bin/*
	cd usr && ln -s lib lib64
#===============================================
endif
