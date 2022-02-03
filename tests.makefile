.PHONY: all clean_temp build compressors storage prime_hints_generator hints modular lbases lthreads sep_fermat special preferred_rules
default: all

EXEC= bin/FIRE6
EXECp= ${EXEC}p
TESTSDIR = ./tests

all: | clean_temp build boxs port compressors storage hints modular lbases lthreads sep_fermat special preferred_rules

build:
	make -f Makefile
	@echo ""

clean_temp:
	rm -rf ${TESTSDIR}/db/*
	rm -rf ${TESTSDIR}/outputs/*
	rm -rf ${TESTSDIR}/hints/*
	rm -rf ${TESTSDIR}/storage/*
	@echo ""

compressors: | clean_temp basiccomp zlib zstd snappy

basiccomp:
	@echo "*** Testing #memory option"
	@echo "*** Testing different compressors"
	${EXEC} -c ${TESTSDIR}/c_none > /dev/null
	diff ${TESTSDIR}/etalon/box.tables ${TESTSDIR}/outputs/c_none.tables
	${EXEC} -c ${TESTSDIR}/c_lz4 > /dev/null
	diff ${TESTSDIR}/etalon/box.tables ${TESTSDIR}/outputs/c_lz4.tables
	${EXEC} -c ${TESTSDIR}/c_lz4fast > /dev/null
	diff ${TESTSDIR}/etalon/box.tables ${TESTSDIR}/outputs/c_lz4fast.tables
	${EXEC} -c ${TESTSDIR}/c_lz4hc > /dev/null
	diff ${TESTSDIR}/etalon/box.tables ${TESTSDIR}/outputs/c_lz4hc.tables
	@echo ""

ifeq ($(findstring --enable-zlib,$(shell cat previous_options)),--enable-zlib)
zlib:
	${EXEC} -c ${TESTSDIR}/c_zlib > /dev/null
	diff ${TESTSDIR}/etalon/box.tables ${TESTSDIR}/outputs/c_zlib.tables
else
zlib:
	@echo "*** zlib compressor option skipped"
endif

ifeq ($(findstring --enable-zstd,$(shell cat previous_options)),--enable-zstd)
zstd:
	${EXEC} -c ${TESTSDIR}/c_zstd > /dev/null
	diff ${TESTSDIR}/etalon/box.tables ${TESTSDIR}/outputs/c_zstd.tables
else
zstd:
	@echo "*** zstd compressor option skipped"
endif

ifeq ($(findstring --enable-snappy,$(shell cat previous_options)),--enable-snappy)
snappy:
	${EXEC} -c ${TESTSDIR}/c_snappy > /dev/null
	diff ${TESTSDIR}/etalon/box.tables ${TESTSDIR}/outputs/c_snappy.tables
else
snappy:
	@echo "*** snappy compressor option skipped"
endif

storage: clean_temp
	rm -rf ${TESTSDIR}/db/*
	rm -rf ${TESTSDIR}/storage/*
	@echo "*** Testing #bucket, #storage options"
	- timeout -s 9 30 ${EXEC} -c ${TESTSDIR}/db_storage >/dev/null
	rm -rf ${TESTSDIR}/db/*
	@echo "*** Dropping FIRE, starting anew on existing storage"
	${EXEC} -c ${TESTSDIR}/db_storage > /dev/null
	diff ${TESTSDIR}/etalon/db_storage.tables ${TESTSDIR}/outputs/db_storage.tables
	@echo ""

port: clean_temp
	@echo "*** Testing #port option and connection from FLAME to FIRE"
	- timeout 30 bin/FLAME6 -c ${TESTSDIR}/port-slave >/dev/null && echo "FLAME" > ${TESTSDIR}/outputs/test & timeout 20 ${EXEC} -c ${TESTSDIR}/port-master >/dev/null
	sleep 5
	diff ${TESTSDIR}/etalon/box.tables ${TESTSDIR}/outputs/port.tables
	echo "FLAME" | diff ${TESTSDIR}/outputs/test -
	@echo ""

prime_hint_generator: clean_temp
	@echo "*** Testing #prime, #small options"
	@echo "*** Generating hints"
	${EXECp} -c ${TESTSDIR}/prime_hint > /dev/null
	diff ${TESTSDIR}/etalon/prime.tables ${TESTSDIR}/outputs/prime.tables
	@echo ""

hints: | clean_temp prime_hint_generator
	@echo "*** Testing #hint option"
	@echo "*** Reading from hint files"
	${EXEC} -c ${TESTSDIR}/hinted > /dev/null
	diff ${TESTSDIR}/etalon/hinted_true.tables ${TESTSDIR}/outputs/hinted.tables
	@echo ""

boxs: clean_temp
	@echo "*** Basic test with global symmetries"
	${EXEC} -c ${TESTSDIR}/boxs > /dev/null
	diff ${TESTSDIR}/etalon/boxs.tables ${TESTSDIR}/outputs/boxs.tables
	@echo ""

special: clean_temp
	@echo "*** Testing #allIBP, #pos_pref options"
	${EXEC} -c ${TESTSDIR}/ibp > /dev/null
	diff ${TESTSDIR}/etalon/box.tables ${TESTSDIR}/outputs/special.tables
	@echo ""

lbases: clean_temp
	@echo "*** Testing #lbases, #wrap option"
	${EXEC} -c ${TESTSDIR}/lbases > /dev/null
	diff ${TESTSDIR}/etalon/lbases.tables ${TESTSDIR}/outputs/lbases.tables
	@echo ""

preferred_rules: clean_temp
	@echo "*** Testing #preferred, #rules options"
	${EXEC} -c ${TESTSDIR}/pref_rules > /dev/null
	diff ${TESTSDIR}/etalon/preferred_rules.tables ${TESTSDIR}/outputs/preferred_rules.tables
	@echo ""

mix: | clean_temp prime_hint_generator
	rm -rf ${TESTSDIR}/db/*
	rm -rf ${TESTSDIR}/storage/*
	@echo "*** Testing mixed options"
	- timeout 15 ${EXEC} -c ${TESTSDIR}/mix > /dev/null
	${EXEC} -c ${TESTSDIR}/mix > /dev/null
	diff ${TESTSDIR}/etalon/preferred_rules.tables ${TESTSDIR}/outputs/mix.tables
	@echo ""

ifeq ($(findstring --enable-lthreads,$(shell cat previous_options)),--enable-lthreads)
lthreads: clean_temp
	@echo "*** Testing lthreads"
	${EXEC} -c ${TESTSDIR}/lthreads > /dev/null
	diff ${TESTSDIR}/etalon/db_storage.tables ${TESTSDIR}/outputs/lthreads.tables
	@echo ""
else
lthreads:
	@echo "*** Skipping lthreads test"
endif

sep_fermat: clean_temp
	@echo "*** Testing separate fermat workers"
	${EXEC} -c ${TESTSDIR}/sep_fermat > /dev/null
	diff ${TESTSDIR}/etalon/db_storage.tables ${TESTSDIR}/outputs/sep_fermat.tables
	@echo ""

modular: clean_temp
	@echo "*** Testing modular calculations"
	${EXECp} -variables 100-100-1 -c ${TESTSDIR}/modular > /dev/null
	diff ${TESTSDIR}/etalon/modular.tables ${TESTSDIR}/outputs/modular-100-100-1.tables
	@echo ""
