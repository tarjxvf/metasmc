#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "tsindex.h"

void tsindex_rebuild(struct tsindex *tr)
{
	double *weights;
	int i;

	weights = malloc(sizeof(double) * tr->maxnodes);
	for(i = 0; i < tr->maxnodes; i++){
		struct edge *e;

		e = tr->edges[i + 1];
		if(e)
			weights[i] = e->top->t - e->bot->t;
		else
			weights[i] = 0;
	}
	__bit_build(tr->index, tr->maxnodes, weights);
	free(weights);
	tsindex_clearflag(tr, TSINDEX_DIRTY);
	tsindex_clearflag(tr, TSINDEX_REBUILD);
}

struct tsindex *tsindex_alloc(int maxedges)
{
	struct tsindex *tr;

	tr = malloc(sizeof(struct tsindex));
	tr->flags = 0;
	tr->index = bit_alloc(maxedges);
	tr->maxedges = maxedges;
	tr->maxnodes = tr->nedges = 0;
	tr->edges = malloc(sizeof(struct edge *) * maxedges);
	memset(tr->edges, 0, sizeof(struct edge *) * maxedges);
	list_init(&tr->free_list);
	list_init(&tr->id_list);

	return tr;
}

void tsindex_reset(struct tsindex *tr)
{
	struct list_head *l;
	struct list *queue;

	tr->flags = 0;
	memset(tr->edges, 0, sizeof(struct edge *) * tr->maxedges);
	tr->nedges = tr->maxnodes = 0;
	bit_clear(tr->index);
	memset(tr->index->freq, 0, sizeof(double) * tr->index->n);

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

struct edge *tsindex_search(struct tsindex *tr, double g, double *cum)
{
	int eid;

//	if(tsindex_dirty(tr)){
	if(tsindex_isrebuild(tr)){
		/* Search operation is disabled when dirty bit is set. */
		return NULL;

	}else{
		eid = bit_getindex(tr->index, g);
		*cum = bit_cumfreq(tr->index, eid - 1);

		return tr->edges[eid];
	}
}

void tsindex_add(struct tsindex *tr, struct edge *e)
{
	struct list *queue;
	double diff;
	int id;

#ifdef DEBUG
	fprintf(stderr, "%s: %d: e=%x(%.6f, %.6f, xtid=%d)\n", __func__, __LINE__, e, e->top->t, e->bot->t, e->xtid);
#endif

	diff = e->top->t - e->bot->t;
	queue = &tr->free_list;

	if(queue->front == NULL){
		/* Allocate a new node in binary indexed tree. */
		if(tsindex_isrebuild(tr)){
//			tsindex_setflag(tr, TSINDEX_DIRTY);
			id = tr->maxnodes + 1;

		}else{
			id = bit_append(tr->index, diff);
		}
		tr->maxnodes++;
		if(tr->maxnodes >= tr->maxedges){
			tr->edges = realloc(tr->edges, sizeof(struct edge *) * (tr->maxedges + 1000));
			memset(tr->edges + tr->maxedges, 0, sizeof(struct edge *) * 1000);
			tr->maxedges += 1000;
		}

	}else{
		struct list_head *l;
		int *ptr;

		/* Get a existing empty node in binary indexed tree. */
		l = queue->front;
		ptr = (int *)GET_OBJ(l);
		id = *ptr;
		__list_remove(queue, l);
		__list_append(&tr->id_list, l);
//		free(l);

		if(!tsindex_isrebuild(tr)){
//			tsindex_setflag(tr, TSINDEX_DIRTY);
//
//		}else{
			bit_update(tr->index, id, diff);
		}
	}

	e->xtid = id;
	tr->edges[id] = e;
	tr->nedges++;
}

void tsindex_update(struct tsindex *tr, struct edge *e, double diff)
{
	if(!tsindex_isrebuild(tr))
//		tsindex_setflag(tr, TSINDEX_DIRTY);
//	else
		bit_update(tr->index, e->xtid, diff);
}

/* Clear a node in binary indexed tree. */
void tsindex_clear(struct tsindex *tr, struct edge *e)
{
	struct list_head *l;
	int *pidx, id;
	double diff;

	id = e->xtid;
	if(!tsindex_isrebuild(tr)){
//		tsindex_setflag(tr, TSINDEX_DIRTY);
//
//	}else{
		diff = bit_getvalue(tr->index, id);
		bit_update(tr->index, id, -diff);
	}

	tr->edges[id] = NULL;
	e->xtid = 0;

	/* Add freed id to the queue. */
	if(tr->id_list.front == NULL)
		l = malloc(sizeof(struct list_head) + sizeof(int));
	else
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

