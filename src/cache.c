#include <stdlib.h>

#include "global.h"
#include "cache.h"

int n_resize = 0;

void cache_resize(struct cache *nc, int add)
{
	int *pid, i, new_size;
	struct intlist *l;
	struct ptrlist *cl;
	char *obj, *objs;

	n_resize++;

	new_size = nc->cache_size + add;
	objs = malloc((nc->obj_size + sizeof(int)) * add);
	cl = malloc(sizeof(struct ptrlist));
	cl->ptr = objs;
	__list_append(&nc->chunk_list, GET_LIST(cl));

	nc->objs = realloc(nc->objs, sizeof(void *) * new_size);
	obj = objs;
	for(i = nc->cache_size; i < new_size; i++){
		nc->objs[i] = obj;
		if(nc->obj_init)
			nc->obj_init(obj, nc->data);
		pid = (int *)(obj + nc->obj_size);
		*pid = i;
		obj += nc->obj_size + sizeof(int);
	}

	for(i = nc->cache_size; i < new_size; i++){
		l = malloc(sizeof(struct intlist));
		__list_append(&nc->id_list, GET_LIST(l));
	}

	nc->cache_size = new_size;
}

void *cache_alloc(struct cache *nc)
{
	struct list *queue;
	struct intlist *l;
	int *ptr, id;

	queue = &nc->free_list;

	if(queue->front){
		/* Get a existing empty node in binary indexed tree. */
		l = (struct intlist *)__list_pop(queue);
		id = l->id;
		__list_append(&nc->id_list, GET_LIST(l));
		return nc->objs[id];

	}else{
		/* Allocate a new node in binary indexed tree. */
		if(nc->maxnodes >= nc->cache_size)
			cache_resize(nc, 1000);
		id = nc->maxnodes++;
		return nc->objs[id];
	}
}

void cache_free(struct cache *nc, void *obj)
{
	struct intlist *l;
	int id;

	id = *(int *)((char *)obj + nc->obj_size);

	/* Add freed id to the queue. */
	l = (struct intlist *)__list_pop(&nc->id_list);
	l->id = id;
	__list_append(&nc->free_list, GET_LIST(l));
}

/* Force all elements to be released.
 * Programmer is reponsible for ensuring that all elements will not be accessed. */
void cache_clear(struct cache *nc)
{
	struct list_head *l;
	struct list *queue;

	queue = &nc->free_list;
	while(queue->front){
		l = __list_pop(queue);
		__list_append(&nc->id_list, l);
	}
	nc->maxnodes = 0;
}

struct cache *cache_create(size_t obj_size, int cache_size, void (*obj_init)(void *obj, void *data), void *data)
{
	struct cache *nc;

	nc = malloc(sizeof(struct cache));

	nc->obj_size = obj_size;
	nc->cache_size = nc->maxnodes = 0;
	nc->objs = NULL;
	if(data)
		nc->data = data;
	else
		nc->data = NULL;

	if(obj_init)
		nc->obj_init = obj_init;
	else
		nc->obj_init = NULL;

	list_init(&nc->chunk_list);
	list_init(&nc->free_list);
	list_init(&nc->id_list);

	cache_resize(nc, cache_size);

	return nc;
}

void cache_destroy(struct cache *nc)
{
	struct list_head *l;
	struct ptrlist *pl;

	cache_clear(nc);
	l = nc->id_list.front;
	while(nc->id_list.front){
		l = __list_pop(&nc->id_list);
		free(l);
	}

	while(nc->chunk_list.front){
		pl = (struct ptrlist *)__list_pop(&nc->chunk_list);
		free(pl->ptr);
		free(pl);
	}

	free(nc->objs);
	free(nc);
}

