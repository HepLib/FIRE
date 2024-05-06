# GMP
GMP:=gmp-6.3.0
GMPx:=tar.gz
# MPFR
MPFR:=mpfr-4.2.1
MPFRx:=tar.gz
# Flint
FLINT:=flint-3.1.3
FLINTx:=tar.gz
# XMalloc
XMALLOC:=jemalloc-5.3.0 # mimalloc-2.1.4 | gperftools-2.10
XMALLOCx:=tar.bz2

ifeq ($(strip $(MALLOC)),)
MALLOC:=-ljemalloc # -lmimalloc | -ltcmalloc
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
cp -r src M Q F P Makefile usr usr-src examples FIREs/; \
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
ifeq ($(UNAME_M),arm64) # Apple Silicon Chip
#===============================================
dep:
	cd usr-src; \
tar xf $(strip $(FLINT)).$(FLINTx); \
tar xf $(strip $(XMALLOC)).$(XMALLOCx);

	cd usr-src/$(FLINT); \
./bootstrap.sh; \
./configure --disable-static --with-gmp=/opt/homebrew --with-mpfr=/opt/homebrew --prefix=$(USR_DIR) CFLAGS="-O3"; \
$(MAKE) install

	cd usr-src/$(XMALLOC); \
./configure --prefix=$(USR_DIR); \
$(MAKE) install

cleandepall: cleandep
	rm -rf usr/include/*
	rm -rf usr/lib/*
	rm -rf usr/share/*
	rm -rf usr/lib64
	rm -rf usr/bin/*
	cd usr && ln -s lib lib64
	
cleandep:
	rm -rf usr-src/$(FLINT)
	rm -rf usr-src/$(XMALLOC)

#===============================================
else
#===============================================
dep:
	cd usr-src; \
tar xf $(strip $(GMP)).$(GMPx); \
tar xf $(strip $(MPFR)).$(MPFRx); \
tar xf $(strip $(FLINT)).$(FLINTx); \
tar xf $(strip $(XMALLOC)).$(XMALLOCx);

	cd usr-src/$(GMP);\
./configure --enable-cxx --prefix=$(USR_DIR); \
$(MAKE) install

	cd usr-src/$(MPFR); \
./configure --with-gmp=$(USR_DIR) --prefix=$(USR_DIR); \
$(MAKE) install

	cd usr-src/$(FLINT); \
./bootstrap.sh; \
./configure --disable-static --enable-avx2 --with-gmp=$(USR_DIR) --with-mpfr=$(USR_DIR) --prefix=$(USR_DIR) CFLAGS="-O3"; \
$(MAKE) install
	
	cd usr-src/$(XMALLOC); \
./configure --prefix=$(USR_DIR); \
$(MAKE) install

cleandep:
	rm -rf usr-src/$(GMP)
	rm -rf usr-src/$(MPFR)
	rm -rf usr-src/$(FLINT)
	rm -rf usr-src/$(XMALLOC)

cleandepall: cleandep
	rm -rf usr/include/*
	rm -rf usr/lib/*
	rm -rf usr/share/*
	rm -rf usr/lib64
	rm -rf usr/bin/*
	cd usr && ln -s lib lib64
#===============================================
endif
