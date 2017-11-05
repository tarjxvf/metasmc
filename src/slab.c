#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "slab.h"

#define PAGE_SIZE	1048576

#define SLAB_STATE_EMPTY	0
#define SLAB_STATE_FULL		1
#define SLAB_STATE_PARTIAL	2
#define SLAB_MAXSTATE		3

#define BUFCTL_STATE_NEW	0
#define BUFCTL_STATE_USED	1
#define BUFCTL_STATE_FREE	2

struct mem_bufctl;

struct bufctl_large {
	union bufctl_small *next;
	struct mem_slab *back;	/* Back pointer to parent slab. */

	char *buf;
};

union bufctl_small {
	union bufctl_small *next;
	struct mem_slab *back;	/* Back pointer to parent slab. */
};

/* Control block of slab. */
struct mem_slab {
	struct mem_slab *prev;
	struct mem_slab *next;

	struct mem_cache *back;	/* Pointer to parent cache. */

	struct mem_bufctl *free_list;

	int freecount;
	int state;

	char *buf;
};

/* Control block of cache */
struct mem_cache {
	size_t slab_size;
	size_t obj_size;
	int obj_per_slab;

#ifdef SLAB_STATS
	int objcount;
	int slabcount;
#endif

	struct mem_slab *slab_lists[SLAB_MAXSTATE];

	/* Function for specific layout. There are two types of layouts. */
	struct mem_slab *(*mem_cache_slab_create)(struct mem_cache *);
	void (*mem_cache_slab_destroy)(struct mem_cache *, struct mem_slab *);

	union bufctl_small *(*mem_cache_findbufctl)(struct mem_cache *, char *);
	void *(*mem_slab_alloc)(struct mem_slab *);
	void (*mem_slab_free)(struct mem_bufctl *);

	void (*constructor)(char *, size_t);
	void (*destructor)(char *, size_t);
};

struct mem_slab *mem_cache_slab_create_small(struct mem_cache *);
void mem_cache_slab_destroy_small(struct mem_cache *, struct mem_slab *);
//__inline__ union bufctl_small *mem_cache_findbufctl_small(struct mem_cache *, char *);
union bufctl_small *mem_cache_findbufctl_small(struct mem_cache *, char *);
void *mem_slab_alloc_small(struct mem_slab *);
void mem_slab_free_small(union bufctl_small *);
void mem_slab_list_remove(struct mem_slab **, struct mem_slab *);
void mem_slab_list_add(struct mem_slab **, struct mem_slab *);
//__inline__ void bufctl_small_list_add(union bufctl_small **, union bufctl_small *);
void bufctl_small_list_add(union bufctl_small **, union bufctl_small *);

struct mem_cache *mem_cache_create(size_t obj_size, void (*constructor)(char *, size_t), void (*destructor)(char *, size_t))
{
	struct mem_cache *new_cache;
	int factor, obj_per_slab;
	size_t unit_size;

	new_cache = malloc(sizeof(struct mem_cache));/* Allocate a new cache */
	memset(new_cache, 0, sizeof(struct mem_cache));

	/* Calculate number of objects per slab,
	 * as well as slab size. */
	unit_size = obj_size + sizeof(struct mem_slab);
	factor = PAGE_SIZE / unit_size;

	/* Calculate slab size. */
	if(factor > 8){
		new_cache->slab_size = PAGE_SIZE;

	}else{
		new_cache->slab_size = 0;

		while(new_cache->slab_size < unit_size * 8)
			new_cache->slab_size += PAGE_SIZE;
	}

	new_cache->obj_per_slab = (new_cache->slab_size - sizeof(struct mem_slab)) /
					(obj_size + sizeof(union bufctl_small));

	/* Set remaining properties of new cache */
	new_cache->obj_size = obj_size;
	new_cache->constructor = constructor;
	new_cache->destructor = destructor;

	/* Use small object layout. */
	new_cache->mem_cache_slab_create = mem_cache_slab_create_small;
	new_cache->mem_cache_slab_destroy = mem_cache_slab_destroy_small;

	new_cache->mem_cache_findbufctl = mem_cache_findbufctl_small;
	new_cache->mem_slab_alloc = mem_slab_alloc_small;
	new_cache->mem_slab_free = mem_slab_free_small;

	/* Allocate new slab. */
	new_cache->slab_lists[SLAB_STATE_EMPTY] = new_cache->mem_cache_slab_create(new_cache);

#ifdef SLAB_STATS
	new_cache->slabcount = 1;
	new_cache->objcount = 0;
#endif
	return new_cache;
}

void *mem_cache_alloc(struct mem_cache *cache)
{
	struct mem_slab *slabctl;
	union bufctl_small *bufctl;
	void *obj;

	/* Find slab that have free object. We start from list of partially full slabs,
	 * Because we tend to reserve empty slabs so that they can be freed if necessary.
	 * If all slabs are full, allocate a new one. */
	if(cache->slab_lists[SLAB_STATE_PARTIAL]){
		slabctl = cache->slab_lists[SLAB_STATE_PARTIAL];

	}else if(cache->slab_lists[SLAB_STATE_EMPTY]){
		slabctl = cache->slab_lists[SLAB_STATE_EMPTY];

	}else{
		/* Get a new slab. */
		slabctl = cache->mem_cache_slab_create(cache);
#ifdef SLAB_STATS
		cache->slabcount++;
#endif
	}

	obj = cache->mem_slab_alloc(slabctl);
#ifdef SLAB_STATS
	cache->objcount++;
#endif

	/* Check if this slab is full. */
	if(slabctl->freecount <= 0){
		/* Move to full list. */
		mem_slab_list_remove(&cache->slab_lists[slabctl->state], slabctl);
		slabctl->state = SLAB_STATE_FULL;
		mem_slab_list_add(&cache->slab_lists[SLAB_STATE_FULL], slabctl);

	}else if(slabctl->state == SLAB_STATE_EMPTY){
		/* If original state is empty, move it to partial list. */
		mem_slab_list_remove(&cache->slab_lists[SLAB_STATE_EMPTY], slabctl);
		slabctl->state = SLAB_STATE_PARTIAL;
		mem_slab_list_add(&cache->slab_lists[SLAB_STATE_PARTIAL], slabctl);
	}

	return obj;
}

void mem_cache_free(struct mem_cache *cache, char *obj)
{
	union bufctl_small *bufctl;
	struct mem_slab *slabctl;

	/* Find the slab containing this object. */
	bufctl = cache->mem_cache_findbufctl(cache, obj);
	slabctl = bufctl->back;
	cache->mem_slab_free(bufctl);
#ifdef SLAB_STATS
	cache->objcount--;
#endif

	/* Check slab state. */
	if(slabctl->freecount >= cache->obj_per_slab){	/* Slab becomes empty */
		/* Insert into empty list. */
		mem_slab_list_remove(&cache->slab_lists[slabctl->state], slabctl);
		slabctl->state = SLAB_STATE_EMPTY;
		mem_slab_list_add(&cache->slab_lists[SLAB_STATE_EMPTY], slabctl);

	}else if(slabctl->state == SLAB_STATE_FULL){	/* Slab is originally full and become partial now. */
		/* Move to partial list. */
		slabctl->state = SLAB_STATE_PARTIAL;
		mem_slab_list_remove(&cache->slab_lists[SLAB_STATE_FULL], slabctl);
		mem_slab_list_add(&cache->slab_lists[SLAB_STATE_PARTIAL], slabctl);
	}
}

/* Destroy a cache */
void mem_cache_destroy(struct mem_cache *cache)
{
	struct mem_slab *slabctl;
	int i;

	for(i = 0; i < SLAB_MAXSTATE; i++){
		slabctl = cache->slab_lists[i];

		while(cache->slab_lists[i]){
			slabctl = cache->slab_lists[i];
			mem_slab_list_remove(&cache->slab_lists[i], slabctl);
			cache->mem_cache_slab_destroy(cache, slabctl);
		}
	}

	/* Release cache itself. */
	free(cache);
}

/****************************************** Functions for layout of small objects. ************************************************/
struct mem_slab *mem_cache_slab_create_small(struct mem_cache *cache)
{
	char *new_slab;
	char *buf;
	struct mem_slab *slabctl;
	union bufctl_small **bufctl;
	size_t unit_size, step_size;
	int i;

	unit_size = cache->obj_size + sizeof(union bufctl_small);

	new_slab = malloc(cache->slab_size);
	memset(new_slab, 0, cache->slab_size);
	slabctl = (struct mem_slab *)(new_slab + cache->slab_size - sizeof(struct mem_slab));	/* Set slab control block. */

	buf = new_slab;

	slabctl->prev = NULL;
	slabctl->next = NULL;
	slabctl->back = cache;
	slabctl->freecount = cache->obj_per_slab;
	slabctl->buf = new_slab;
	slabctl->state = SLAB_STATE_EMPTY;

	bufctl = &slabctl->free_list;

	/* Initialize all objects in new slab. */
	step_size = cache->obj_size + sizeof(union bufctl_small);

	if(cache->constructor){
		for(i = 0; i < cache->obj_per_slab; i++){
			*bufctl = (union bufctl_small *)(buf + cache->obj_size);
			cache->constructor(buf, cache->obj_size);
			buf += step_size;
			bufctl = &(*bufctl)->next;
		}

	}else{
		for(i = 0; i < cache->obj_per_slab; i++){
			*bufctl = (union bufctl_small *)(buf + cache->obj_size);
			buf += step_size;
			bufctl = &(*bufctl)->next;
		}
	}

	/* Set freelist point of last object to NULL. */
	*bufctl = NULL;

	return slabctl;
}

void mem_cache_slab_destroy_small(struct mem_cache *cache, struct mem_slab *slabctl)
{
	char *buf;
	union bufctl_small *bufctl;
	size_t step_size;
	int i;

	step_size = cache->obj_size + sizeof(union bufctl_small);

	buf = slabctl->buf;
	bufctl = (union bufctl_small *)(buf + cache->obj_size);

	/* Destruct all objects. */
	for(i = 0; i < cache->obj_per_slab; i++){
		if(cache->destructor)
			cache->destructor(buf, cache->obj_size);

		buf += step_size;
		bufctl = (union bufctl_small *)(buf + cache->obj_size);
	}

	/* Remove slab from list. */
	mem_slab_list_remove(&cache->slab_lists[slabctl->state], slabctl);

	/* Free slab memory */
	free(slabctl->buf);
//	cache->slabcount--;
}

//__inline__ union bufctl_small *mem_cache_findbufctl_small(struct mem_cache *cache, char *obj)
union bufctl_small *mem_cache_findbufctl_small(struct mem_cache *cache, char *obj)
{
	return (union bufctl_small *)(obj + cache->obj_size);
}

void *mem_slab_alloc_small(struct mem_slab *slabctl)
{
	char *obj;
	union bufctl_small *bufctl;

	/* Find first object in free list */
	bufctl = slabctl->free_list;
	slabctl->free_list = bufctl->next;
	bufctl->back = slabctl;
	slabctl->freecount--;

	return (char *)bufctl - slabctl->back->obj_size;
}

void mem_slab_free_small(union bufctl_small *bufctl)
{
	struct mem_slab *slabctl;

	slabctl = bufctl->back;
	/* Insert into free list. */
	bufctl_small_list_add(&slabctl->free_list, bufctl);
	/* Increate slab freecount. */
	slabctl->freecount++;
}

void mem_slab_list_remove(struct mem_slab **list_head, struct mem_slab *slabctl)
{
	if(slabctl->prev)
		slabctl->prev->next = slabctl->next;

	if(slabctl->next)
		slabctl->next->prev = slabctl->prev;

	if(*list_head == slabctl)
		*list_head = slabctl->next;
}

void mem_slab_list_add(struct mem_slab **list, struct mem_slab *slabctl)
{
	slabctl->next = *list;
	slabctl->prev = NULL;

	if(*list)
		(*list)->prev = slabctl;

	*list = slabctl;
}

//__inline__ void bufctl_small_list_add(union bufctl_small **list, union bufctl_small *bufctl)
void bufctl_small_list_add(union bufctl_small **list, union bufctl_small *bufctl)
{
	bufctl->next = *list;
	*list = bufctl;
}

void mem_cache_release(struct mem_cache *cache)
{
	struct mem_slab *slabctl;

	while(cache->slab_lists[SLAB_STATE_EMPTY]){
		slabctl = cache->slab_lists[SLAB_STATE_EMPTY];
		mem_slab_list_remove(&cache->slab_lists[SLAB_STATE_EMPTY], slabctl);
		cache->mem_cache_slab_destroy(cache, slabctl);
	}
}

/* Print statistics about memory usage. */
void mem_cache_stats(struct mem_cache *cache)
{
	printf("Number of slabs: %d\n", cache->slabcount);
	printf("Slab size: %d\n", cache->slab_size);
	printf("Number of objects used: %d\n", cache->objcount);
	printf("Object size: %d\n", cache->obj_size);
}
