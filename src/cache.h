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
	struct list chunk_list;
	struct list free_list;	// List of free index before maxnodes
	struct list id_list;

	void *data;
	void (*obj_init)(void *obj, void *data);
};

void cache_resize(struct cache *nc, int add);
void *cache_alloc(struct cache *nc);
void cache_free(struct cache *nc, void *obj);

void cache_clear(struct cache *nc);
struct cache *cache_create(size_t obj_size, int cache_size, void (*obj_init)(void *obj, void *data), void *data);
void cache_destroy(struct cache *nc);

extern int n_resize;

#endif
