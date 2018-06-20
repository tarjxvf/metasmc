#include <stdlib.h>

#include "global.h"
#include "cache.h"

void cache_resize(struct cache *nc, int add)
{
	int *pid, i, new_size;
	struct list_head *l, *cl;
	char *obj, *objs;

	new_size = nc->cache_size + add;
	objs = malloc((nc->obj_size + sizeof(int)) * add);
	cl = malloc(sizeof(char *) + sizeof(struct list_head));
	*( (char **) GET_OBJ(cl)) = objs;
	__list_append(&nc->chunk_list, cl);

	nc->objs = realloc(nc->objs, sizeof(void *) * new_size);
	obj = objs;
	for(i = nc->cache_size; i < new_size; i++){
//		nc->objs[i] = obj = malloc(nc->obj_size + sizeof(int));
		nc->objs[i] = obj;
		pid = (int *)(obj + nc->obj_size);
		*pid = i;
		obj += nc->obj_size + sizeof(int);
	}

	for(i = nc->cache_size; i < new_size; i++){
		l = malloc(sizeof(struct list_head) + sizeof(int));
		__list_append(&nc->id_list, l);
	}

	nc->cache_size = new_size;
}

void *cache_alloc(struct cache *nc)
{
	struct list *queue;
	int id;

	queue = &nc->free_list;

	if(queue->front == NULL){
		/* Allocate a new node in binary indexed tree. */
		if(nc->maxnodes >= nc->cache_size){
			cache_resize(nc, 1000);
		}
		id = nc->maxnodes++;

	}else{
		struct list_head *l;
		int *ptr;

		/* Get a existing empty node in binary indexed tree. */
		l = queue->front;
		ptr = (int *)GET_OBJ(l);
		id = *ptr;
		__list_remove(queue, l);
		__list_append(&nc->id_list, l);
//		free(l);
	}

	return nc->objs[id];
}

void cache_free(struct cache *nc, void *obj)
{
	struct list_head *l;
	int id, *pidx;

	id = *(int *)((char *)obj + nc->obj_size);

	/* Add freed id to the queue. */
//	if(nc->id_list.front == NULL)
//		l = malloc(sizeof(struct list_head) + sizeof(int));
//	else
		l = __list_pop(&nc->id_list);
	pidx = (int *)GET_OBJ(l);
	*pidx = id;
	__list_append(&nc->free_list, l);
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

struct cache *cache_create(size_t obj_size, int cache_size)
{
	struct list_head *l;
	struct cache *nc;
	int i, *pid;
	char *obj;

	nc = malloc(sizeof(struct cache));

	nc->obj_size = obj_size;
//	nc->cache_size = cache_size;
	nc->cache_size = nc->maxnodes = 0;
	nc->objs = NULL;

	list_init(&nc->chunk_list);
	list_init(&nc->free_list);
	list_init(&nc->id_list);

	cache_resize(nc, cache_size);

/*	nc->objs = malloc(sizeof(void *) * cache_size);
	for(i = 0; i < cache_size; i++){
		obj = (char *)malloc(obj_size + sizeof(int));
		nc->objs[i] = (void *)obj;
		pid = (int *)(obj + obj_size);
		*pid = i;
	}

	for(i = 0; i < nc->cache_size; i++){
		l = malloc(sizeof(struct list_head) + sizeof(int));
		__list_append(&nc->id_list, l);
	}*/

	return nc;
}

void cache_destroy(struct cache *nc)
{
	struct list_head *l;
	char *obj;
	int i;

	cache_clear(nc);
	l = nc->id_list.front;
	while(nc->id_list.front){
		l = __list_pop(&nc->id_list);
		free(l);
	}

//	for(i = 0; i < nc->cache_size; i++)
//		free(nc->objs[i]);

	while(nc->chunk_list.front){
		l = __list_pop(&nc->chunk_list);
		obj = *((char **)GET_OBJ(l));
		free(obj);
		free(l);
	}

	free(nc->objs);
	free(nc);
}

