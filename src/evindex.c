#include <stdlib.h>
#include <math.h>

#include "evindex.h"

/***** LIBAVL allocator of population-wise event index. *****/
int evindex_compar(struct event *a, struct event *b)
{
	double diff;

	diff = a->t - b->t;
	if(diff > 0){
		return 1;

	}else if(diff < 0){
		return -1;

	}else{
		return a->type - b->type;
	}
}

void evindex_seq_off(struct evindex *evidx)
{
	struct rb_node **nodes;
	int nnodes, i, j;

	nnodes = evidx->idx->ls.n;
	nodes = malloc(sizeof(struct rb_node *) * nnodes);
	rbindex_clearflag(evidx->idx, RBINDEX_SEQUENTIAL);
	rbindex_rebuild_tree(evidx->idx, nodes);

	for(i = nnodes - 1; i >= 0; i--){
		struct event *ev, *left, *right;

		ev = nodes[i]->rb_data;
		for(j = 0; j < evidx->npop_all; j++)
			ev->sumdn[j] = ev->dn[j];

		// Link nodes
		if(i * 2 + 1 < nnodes){
			left = nodes[i * 2 + 1]->rb_data;
			for(j = 0; j < evidx->npop_all; j++)
				ev->sumdn[j] += left->sumdn[j];
		}

		if(i * 2 + 2 < nnodes){
			right = nodes[i * 2 + 2]->rb_data;
			for(j = 0; j < evidx->npop_all; j++)
				ev->sumdn[j] += right->sumdn[j];
		}
	}

	free(nodes);
}

void evindex_propagate(struct rb_traverser *tr, int npop, int *dn)
{
	struct event *ev;
	int i, j;

	// Backtrack and update summary statistics
	for(i = tr->rb_height - 1; i >= 0; i++){
		ev = tr->rb_stack[i]->rb_data;
		for(j = 0; j < npop; j++)
			ev->sumdn[j] += dn[j];
	}
}

void evindex_delete(struct evindex *evidx, struct event *ev)
{
	if(rbindex_isseq(evidx->idx)){
		rbindex_s_delete(evidx->idx, (void *)ev);

	}else{
//		if(ev->type == EVENT_JOIN || ev->type == EVENT_SPLT){
//		}else{
			evindex_rb_delete(evidx, ev);
//		}
	}
}

void evindex_insert(struct evindex *evidx, struct event *ev)
{
	if(rbindex_isseq(evidx->idx)){
		evindex_s_insert(evidx, ev);

	}else{
/*		rb_traverser tr;
		struct event *ev;
		int i, j;

		rb_t_insert(&tr, evidx->idx, ev);
		trav_refresh(&tr);
		evindex_propgate(&tr, evidx->npop_all, ev->dn);*/

		evindex_rb_insert(evidx, ev);
	}
}

void evindex_destroy(struct genealogy *G, struct evindex *evidx)
{
	struct event *ev;

	ev = (struct event *)GET_OBJ(evidx->idx->lsentinel);
	free(GET_LIST(ev));
	ev = (struct event *)GET_OBJ(evidx->idx->rsentinel);
	free(GET_LIST(ev));

	rbindex_destroy(evidx->idx);
	free(evidx->dn);
	free(evidx);
}

// Sequential seek for interval endpoint.
void evindex_s_seek(struct evindex *evidx, double t)
{
	if(rbindex_isseq(evidx->idx)){
		seq_traverser lfwd, lbwd;
		struct event *ebwd, *efwd;

		lbwd = evidx->idx->cur_s;
		ebwd = evindex_cur(lbwd);
		while(ebwd->t >= t) ebwd = evindex_prev(&lbwd);
		lfwd = lbwd;
		efwd = evindex_cur(lfwd);
		while(efwd->t < t) efwd = evindex_next(&lfwd);
		evidx->idx->cur_s = lfwd;
	}
}

struct evindex *evindex_create(struct genealogy *G, struct config *cfg)
{
	struct evindex *evidx;
	struct event *ev;
	int npop;

	evidx = malloc(sizeof(struct evindex));
	evidx->idx = rbindex_create(evindex_compar, NULL, &rbindex_allocator);
	evidx->npop_all = npop = cfg->npop_all;
	evidx->dn = malloc(sizeof(int) * npop);

	// Set up sentinel
	ev = alloc_event(G->cfg, EVENT_SAMP, 0);
	evidx->idx->lsentinel = GET_LIST(ev);

	list_append(&evidx->idx->ls, ev);

	ev = alloc_event(G->cfg, EVENT_GSIZ, INFINITY);
	((struct gsiz_event *)ev)->size = 1;
	evidx->idx->rsentinel = GET_LIST(ev);

	list_append(&evidx->idx->ls, ev);

	return evidx;
}

