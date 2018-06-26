#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "tsindex.h"

#define MAXNSEC 1000000000

unsigned long long t_ts_rebuild = 0;
unsigned long long n_ts_resize = 0;

void tsindex_resize(struct tsindex *tr, int add)
{
	struct list_head *l;
	int i;

	n_ts_resize++;

	tr->edges = realloc(tr->edges, sizeof(struct edge *) * (tr->maxedges + add));
	memset(tr->edges + tr->maxedges, 0, sizeof(struct edge *) * add);

	tr->weights = realloc(tr->weights, sizeof(double) * (tr->maxedges + add));
	memset(tr->weights + tr->maxedges, 0, sizeof(double) * add);

	for(i = 0; i < add; i++){
		l = malloc(sizeof(struct list_head) + sizeof(int));
		__list_append(&tr->id_list, l);
	}

	tr->maxedges += add;
}

void tsindex_rebuild(struct tsindex *tr)
{
	double *weights;
	int i;

	struct timespec beg, end;
	int nsec;

//	clock_gettime(CLOCK_MONOTONIC, &beg);

//	weights = malloc(sizeof(double) * tr->maxnodes);
	weights = tr->weights + 1;
	for(i = 0; i < tr->maxnodes; i++){
		struct edge *e;

		e = tr->edges[i + 1];
		if(e)
			weights[i] = e->top->t - e->bot->t;
		else
			weights[i] = 0;
	}
	__bit_build(tr->index, tr->maxnodes, weights);
//	free(weights);
	tsindex_clearflag(tr, TSINDEX_DIRTY);
	tsindex_clearflag(tr, TSINDEX_REBUILD);

//	clock_gettime(CLOCK_MONOTONIC, &end);
//	nsec = (end.tv_sec - beg.tv_sec) * MAXNSEC + (end.tv_nsec - beg.tv_nsec);
//	t_ts_rebuild += nsec;
}

struct tsindex *tsindex_alloc(int maxedges)
{
	struct list_head *l;
	struct tsindex *tr;
	int i;

	tr = malloc(sizeof(struct tsindex));
	tr->flags = 0;
	tr->index = bit_alloc(maxedges);
//	tr->maxedges = 0;
	tr->maxedges = 1;
	tr->maxnodes = tr->nedges = 0;
	tr->edges = malloc(sizeof(struct edge *));
	tr->edges[0] = NULL;
	tr->weights = malloc(sizeof(double));
	tr->weights[0] = 0;
	list_init(&tr->free_list);
	list_init(&tr->id_list);

	tsindex_resize(tr, maxedges);

	return tr;
}

void tsindex_reset(struct tsindex *tr)
{
	struct list_head *l;
	struct list *queue;

	tr->flags = 0;
//	memset(tr->edges, 0, sizeof(struct edge *) * tr->maxedges);
	tr->nedges = tr->maxnodes = 0;
	bit_clear(tr->index);
//	memset(tr->index->freq, 0, sizeof(double) * tr->index->n);

	/* Clear free list. */
	queue = &tr->free_list;
	l = queue->front;
	while(l){
		struct list_head *tmp;

		tmp = l->next;
		__list_remove(queue, l);
		__list_append(&tr->id_list, l);
//		free(l);
		l = tmp;
	}
}

// Do not update the tree.
void tsindex_add_rebuild(struct tsindex *tr, struct edge *e)
{
	struct list_head *l;
	int id, *ptr, i;

	if(tr->free_list.front == NULL){
		/* Allocate a new node in binary indexed tree. */
		id = tr->maxnodes + 1;
		tr->maxnodes++;
		if(tr->maxnodes >= tr->maxedges)
			tsindex_resize(tr, tr->maxedges);

	}else{
		/* Get a existing empty node in binary indexed tree. */
//		l = __list_pop(queue);
		l = tr->free_list.front;
		ptr = (int *)GET_OBJ(l);
		id = *ptr;
		__list_remove(&tr->free_list, l);
		__list_append(&tr->id_list, l);
	}

	e->xtid = id;
	tr->edges[id] = e;
	tr->nedges++;
}

/* Clear a node without updating the tree. */
void tsindex_clear_rebuild(struct tsindex *tr, struct edge *e)
{
	struct list_head *l;
	int *pidx, id;

	id = e->xtid;
	tr->edges[id] = NULL;
	e->xtid = 0;

	/* Add freed id to the queue. */
	l = __list_pop(&tr->id_list);
	pidx = (int *)GET_OBJ(l);
	*pidx = id;
	__list_append(&tr->free_list, l);
	tr->nedges--;
}

void tsindex_add(struct tsindex *tr, struct edge *e)
{
	struct list_head *l;
	double diff;
	int id, *ptr, i;


#ifdef DEBUG
	fprintf(stderr, "%s: %d: e=%x(%.6f, %.6f, xtid=%d)\n", __func__, __LINE__, e, e->top->t, e->bot->t, e->xtid);
#endif

	diff = e->top->t - e->bot->t;

	if(tr->free_list.front == NULL){
		/* Allocate a new node in binary indexed tree. */
		if(tsindex_isrebuild(tr)){
			id = tr->maxnodes + 1;

		}else{
			id = bit_append(tr->index, diff);
		}
		tr->maxnodes++;
		if(tr->maxnodes >= tr->maxedges)
			tsindex_resize(tr, tr->maxedges);

	}else{
		/* Get a existing empty node in binary indexed tree. */
//		l = __list_pop(queue);
		l = tr->free_list.front;
		ptr = (int *)GET_OBJ(l);
		id = *ptr;
		__list_remove(&tr->free_list, l);
		__list_append(&tr->id_list, l);

		if(!tsindex_isrebuild(tr)){
			bit_update(tr->index, id, diff);
		}
	}

	e->xtid = id;
	tr->edges[id] = e;
	tr->nedges++;
}

/* Clear a node in binary indexed tree. */
void tsindex_clear(struct tsindex *tr, struct edge *e)
{
	struct list_head *l;
	int *pidx, id;
	double diff;

	id = e->xtid;
#ifdef DEBUG
	fprintf(stderr, "%s: %d: id=%d\n", __func__, __LINE__, id);
#endif
	if(!tsindex_isrebuild(tr)){
		diff = tr->edges[id]->top->t - tr->edges[id]->bot->t;
		bit_update(tr->index, id, -diff);
	}

	tr->edges[id] = NULL;
	e->xtid = 0;

	/* Add freed id to the queue. */
	l = __list_pop(&tr->id_list);
	pidx = (int *)GET_OBJ(l);
	*pidx = id;
	__list_append(&tr->free_list, l);
	tr->nedges--;
}

void tsindex_free(struct tsindex *tr)
{
	struct list_head *l;

	while(tr->free_list.front){
		l = __list_pop(&tr->free_list);
		free(l);
	}

	while(tr->id_list.front){
		l = __list_pop(&tr->id_list);
		free(l);
	}

	free(tr->edges);
	free(tr->weights);
	bit_free(tr->index);
	free(tr);
}

void tsindex_dump(struct tsindex *tr)
{
	struct list_head *l;
	struct edge *e;
	int i;

	fprintf(stderr, "tr->flags=%d, tr->nedges=%d, tr->maxnodes=%d, tr->maxedges=%d, G->tr_xoveredges=", tr->flags, tr->nedges, tr->maxnodes, tr->maxedges);
	for(i = 1; i <= tr->maxnodes; i++){
		e = tr->edges[i];
		if(e)
			fprintf(stderr, "%x(%.10f, %.10f, xtid=%d), ", e, e->bot->t, e->top->t, e->xtid);
		else
			fprintf(stderr, "NULL, ");
	}
	fprintf(stderr, "\n");

	fprintf(stderr, "free_list=");
	l = tr->free_list.front;
	while(l){
		int *pidx;
		pidx = (int *)GET_OBJ(l);
		fprintf(stderr, "%d->", *pidx);
		l = l->next;
	}
	fprintf(stderr, "\n");

	bit_print(tr->index);
}

