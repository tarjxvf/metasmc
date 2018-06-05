#ifndef EVINDEX_H
#define EVINDEX_H

#include "rbindex.h"
#include "smc.h"
#include "global.h"

struct evindex {
	struct rbindex *idx;
	int npop_all;
	int *dn;
};

static inline struct event *evindex_cur(rb_traverser it)
{
	return (struct event *)rbindex_cur(it);
}

static inline struct event *evindex_prev(rb_traverser *it)
{
	return (struct event *)rbindex_prev(it);
}

static inline struct event *evindex_next(rb_traverser *it)
{
	return (struct event *)rbindex_next(it);
}

static inline struct event *evindex_s_get(struct evindex *evidx)
{
	return (struct event *)GET_OBJ(evidx->idx->cur_s);
}

static inline struct event *evindex_s_forward(struct evindex *evidx)
{
	return (struct event *)evindex_next(&evidx->idx->cur_s);
}

static inline void evindex_s_set(struct evindex *evidx, struct event *e)
{
	rbindex_s_set(evidx->idx, e);
}

void evindex_s_seek(struct evindex *evidx, double t);

static inline void evindex_s_rewind(struct evindex *evidx)
{
	evidx->idx->cur_s = evidx->idx->ls.front->next;
//	evidx->cur_s = evidx->ls.front;
}

static inline int evindex_s_final(struct evindex *evidx)
{
	evidx->idx->cur_s = evidx->idx->rsentinel;
}

static inline int evindex_s_end(struct evindex *evidx)
{
	return evidx->idx->cur_s == NULL;
}

static inline void evindex_seq_on(struct evindex *evidx)
{
	rbindex_seq_on(evidx->idx);
}

void evindex_seq_off(struct evindex *evidx);

static inline struct event *evindex_find(rb_traverser *it, struct evindex *evidx, struct event *key)
{
	return rbindex_find(it, evidx->idx, key);
}

static inline void evindex_s_delete(struct evindex *evidx, struct event *e)
{
	rbindex_s_delete(evidx->idx, (void *)e);
}

static inline void evindex_rb_delete(struct evindex *evidx, struct event *e)
{
	rbindex_rb_delete(evidx->idx, (void *)e);
}

static inline void evindex_delete(struct evindex *evidx, struct event *e)
{
	rbindex_delete(evidx->idx, (void *)e);
}

static inline void evindex_rb_insert(struct evindex *evidx, struct event *e)
{
	rbindex_rb_insert(evidx->idx, (void *)e);
}

static inline void evindex_s_insert(struct evindex *evidx, struct event *e)
{
	rbindex_s_insert(evidx->idx, (void *)e);
}

static inline void evindex_insert(struct evindex *evidx, struct event *e)
{
	if(rbindex_isseq(evidx->idx))
		evindex_s_insert(evidx, e);
	else
		evindex_rb_insert(evidx, e);
}

void evindex_destroy(struct genealogy *G, struct evindex *evidx);
struct evindex *evindex_create(struct genealogy *G, struct config *cfg);

#endif
