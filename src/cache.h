#ifndef CACHE_H
#define CACHE_H

#include "list.h"

struct list;
//#include "global.h"

struct cache {
	int obj_size;
	int cache_size;	// Number of elements in the cacne
	int maxnodes;	// The index of rightmost elements + 1
	void **objs;
	struct list free_list;	// List of free index before maxnodes
	struct list id_list;
};

void cache_resize(struct cache *nc, int add);
void *cache_alloc(struct cache *nc);
void cache_free(struct cache *nc, void *obj);

void cache_clear(struct cache *nc);
struct cache *cache_create(size_t obj_size, int cache_size);
void cache_destroy(struct cache *nc);

#endif
