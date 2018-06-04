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
		return 0;
	}
}

void evindex_seq_off(struct rbindex *evidx)
{
	struct rb_node **nodes;
	int nnodes;

	nnodes = evidx->ls.n;
	nodes = malloc(sizeof(struct rb_node *) * nnodes);
	rbindex_clearflag(evidx, RBINDEX_SEQUENTIAL);
	rbindex_rebuild_tree(evidx, nodes);
	free(nodes);
}

void evindex_destroy(struct genealogy *G, struct rbindex *evidx)
{
	struct event *ev;

	ev = (struct event *)GET_OBJ(evidx->lsentinel);
	free(GET_LIST(ev));
	ev = (struct event *)GET_OBJ(evidx->rsentinel);
	free(GET_LIST(ev));

	rbindex_destroy(evidx);
}

// Sequential seek for interval endpoint.
void evindex_s_seek(struct rbindex *evidx, double t)
{
	if(rbindex_isseq(evidx)){
		rb_traverser lfwd, lbwd;
		struct event *ebwd, *efwd;

		lbwd = evidx->cur_s;
		ebwd = evindex_cur(lbwd);
		while(ebwd->t >= t) ebwd = evindex_prev(&lbwd);
		lfwd = lbwd;
		efwd = evindex_cur(lfwd);
		while(efwd->t < t) efwd = evindex_next(&lfwd);
		evidx->cur_s = lfwd;
	}
}

struct rbindex *evindex_create(struct genealogy *G, struct config *cfg)
{
	struct rbindex *evidx;
	struct event *ev;

	evidx = rbindex_create(evindex_compar, NULL, &rbindex_allocator);

	// Set up sentinel
	ev = alloc_event(G->cfg, EVENT_SAMP, 0);
	evidx->lsentinel = GET_LIST(ev);

	list_append(&evidx->ls, ev);

	ev = alloc_event(G->cfg, EVENT_GSIZ, INFINITY);
	((struct gsiz_event *)ev)->size = 1;
	evidx->rsentinel = GET_LIST(ev);

	list_append(&evidx->ls, ev);

	return evidx;
}

