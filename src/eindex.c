#include <stdlib.h>
#include <math.h>

#include "eindex.h"

/***** LIBAVL allocator of population-wise edge index. *****/
int eindex_compar(struct edge *a, struct edge *b)
{
/*	double diff, diff2;

	diff = a->top->t - b->top->t;
	if(diff > 0){
		return 1;

	}else if(diff < 0){
		return -1;

	}else{
		diff2 = a->bot->t - b->bot->t;
		if(diff2 < 0)
			return -1;
		else if(diff2 > 0)
			return 1;
		else
			return a->eid - b->eid;
	}*/
}

void eindex_seq_off(struct rbindex *eidx)
{
/*	struct rb_node **nodes;
	int nnodes;

	nnodes = eidx->ls.n;
	nodes = eidx->nc.nodes;
	rbindex_clearflag(eidx, RBINDEX_SEQUENTIAL);
	rbindex_rebuild_tree(eidx);*/
}

void eindex_reset(struct genealogy *G, struct rbindex *eidx)
{
/*	list_init(&eidx->ls);
	list_append(&eidx->ls, GET_OBJ(eidx->lsentinel));
	list_append(&eidx->ls, GET_OBJ(eidx->rsentinel));
	rbindex_rb_clear(eidx);*/
}

// Sequential seek for interval endpoint.
void eindex_s_seek_ttop(struct rbindex *eidx, double ttop)
{
/*	if(rbindex_isseq(eidx)){
		seq_traverser lfwd, lbwd;
		struct edge *ebwd, *efwd;

		lbwd = eidx->cur_s;
		ebwd = eindex_cur(lbwd);
		while(ebwd->top->t >= ttop) ebwd = eindex_prev(&lbwd);
		lfwd = lbwd;
		efwd = eindex_cur(lfwd);
		while(efwd->top->t < ttop) efwd = eindex_next(&lfwd);
		eidx->cur_s = lfwd;
	}*/
}

// Sequential seek for time interval
void eindex_s_seek_t(struct rbindex *eidx, double ttop, double tbot)
{
/*	if(rbindex_isseq(eidx)){
		seq_traverser lfwd, lbwd;
		struct edge *ebwd, *efwd;

		lbwd = eidx->cur_s;
		ebwd = eindex_cur(lbwd);
		while(ebwd->top->t > ttop) ebwd = eindex_prev(&lbwd);
		while(ebwd->top->t == ttop && ebwd->bot->t > tbot) ebwd = eindex_prev(&lbwd);

		lfwd = lbwd;
		efwd = eindex_cur(lfwd);
		while(efwd->top->t < ttop) efwd = eindex_next(&lfwd);
		while(efwd->top->t == ttop && ebwd->bot->t < tbot) efwd = eindex_next(&lfwd);

		eidx->cur_s = lfwd;
	}*/
}

// Sequential seek for full key
void eindex_s_seek(struct rbindex *eidx, double ttop, double tbot, int eid)
{
/*	if(rbindex_isseq(eidx)){
		seq_traverser lfwd, lbwd;
		struct edge *ebwd, *efwd;

		lbwd = eidx->cur_s;
		ebwd = eindex_cur(lbwd);
		while(ebwd->top->t > ttop) ebwd = eindex_prev(&lbwd);
		while(ebwd->top->t == ttop && ebwd->bot->t > tbot) ebwd = eindex_prev(&lbwd);
		while(ebwd->top->t == ttop && ebwd->bot->t == tbot && ebwd->eid > eid) ebwd = eindex_prev(&lbwd);

		lfwd = lbwd;
		efwd = eindex_cur(lfwd);
		while(efwd->top->t < ttop) efwd = eindex_next(&lfwd);
		while(efwd->top->t == ttop && efwd->bot->t < tbot) efwd = eindex_next(&lfwd);
		while(efwd->top->t == ttop && efwd->bot->t == tbot && efwd->eid < eid) efwd = eindex_next(&lfwd);

		eidx->cur_s = lfwd;
	}*/
}

void eindex_destroy(struct genealogy *G, struct rbindex *eidx)
{
	struct edge *e;

	e = (struct edge *)GET_OBJ(eidx->lsentinel);
	free_node(G, e->top);
	free_node(G, e->bot);
	free_edge(G, e);

	e = (struct edge *)GET_OBJ(eidx->rsentinel);
	free_node(G, e->top);
	free_node(G, e->bot);
	free_edge(G, e);

	rbindex_destroy(eidx);
}

struct rbindex *eindex_create(struct genealogy *G, int pop)
{
	struct rbindex *eidx;
	struct node *top, *bot;
	struct edge *e;

	eidx = rbindex_create(eindex_compar, G->cfg->maxfrag * 2);

	// Set up sentinel
	top = alloc_node(G, NODE_DUMMY, pop, 0);
	bot = alloc_node(G, NODE_DUMMY, pop, 0);
	e = alloc_edge(G, top, bot);
	eidx->lsentinel = GET_LIST(e);

	list_append(&eidx->ls, e);

	top = alloc_node(G, NODE_DUMMY, pop, INFINITY);
	bot = alloc_node(G, NODE_DUMMY, pop, INFINITY);
	e = alloc_edge(G, top, bot);
	eidx->rsentinel = GET_LIST(e);

	list_append(&eidx->ls, e);

	return eidx;
}

