
AUTOMAKE_OPTIONS=foreign subdir-objects

bin_PROGRAMS=sampler metasmc metaprof
lib_LTLIBRARIES=libsmc.la
include_HEADERS=util.h list.h global.h mutation.h smc.h cache.h bit.h rb.h rbindex.h tsindex.h evindex.h

libsmc_la_SOURCES=smc.c util.c list.c mutation.c global.c rand.c cache.c bit.c tsindex.c rb.c rbindex.c evindex.c dSFMT.c
libsmc_la_LDFLAGS= -shared -fPIC $(LAPACK_LIBS) -lrt -lm -lgfortran -llapack -llapacke
libsmc_la_CFLAGS=$(LAPACK_CFLAGS) 

sampler_SOURCES=sampler.c
sampler_LDADD=libsmc.la
sampler_LDFLAGS=-lrt -lm

metaprof_SOURCES=metaprof.c
metaprof_LDADD=libsmc.la
metaprof_LDFLAGS=-lrt -lm

metasmc_SOURCES=simulator.c
metasmc_LDADD=libsmc.la
metasmc_LDFLAGS=-lrt -lm
