#!/bin/sh

prefix=/Volumes/DATA/Projects/FIREs/usr
exec_prefix=/Volumes/DATA/Projects/FIREs/usr
libdir=${exec_prefix}/lib

DYLD_INSERT_LIBRARIES=${libdir}/libjemalloc.2.dylib
export DYLD_INSERT_LIBRARIES
exec "$@"
