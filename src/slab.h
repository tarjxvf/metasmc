#ifndef SLAB_H
#define SLAB_H

#include <stdio.h>
#define SLAB_STATS

struct mem_bufctl;
struct mem_slab;
struct mem_cache;

/* User interface to cache */
struct mem_cache *mem_cache_create(size_t, void (*)(char *, size_t), void (*)(char *, size_t));
/* Get a free object from cache. */
void *mem_cache_alloc(struct mem_cache *);
/* Return an object back to cache. */
void mem_cache_free(struct mem_cache *, char *);

/* Release memory of free objects back to system. */
void mem_cache_release(struct mem_cache *);
void mem_cache_destroy(struct mem_cache *);
void mem_cache_stats(struct mem_cache *);

#endif
