
AUTOMAKE_OPTIONS=foreign subdir-objects

bin_PROGRAMS=sampler metasmc abundance
lib_LTLIBRARIES=libsmc.la
include_HEADERS=util.h global.h mutation.h smc.h slab.h

libsmc_la_SOURCES=smc.c util.c mutation.c global.c rand.c slab.c mt19937-64/mt19937-64.c
libsmc_la_LDFLAGS= -shared -fPIC $(LAPACK_LIBS) -lrt -lm -lgfortran -llapack -llapacke
libsmc_la_CFLAGS=$(LAPACK_CFLAGS)

sampler_SOURCES=sampler.c
sampler_LDADD=libsmc.la
sampler_LDFLAGS=-lrt -lm

abundance_SOURCES=abundance.c
abundance_LDADD=libsmc.la
abundance_LDFLAGS=-lrt -lm

metasmc_SOURCES=simulator.c
metasmc_LDADD=libsmc.la
metasmc_LDFLAGS=-lrt -lm
