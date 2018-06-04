#ifndef EVINDEX_H
#define EVINDEX_H

#include "rbindex.h"
#include "smc.h"
#include "global.h"

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

static inline struct event *evindex_s_get(struct rbindex *evidx)
{
	return (struct event *)GET_OBJ(evidx->cur_s);
}

static inline struct event *evindex_s_forward(struct rbindex *evidx)
{
	return (struct event *)evindex_next(&evidx->cur_s);
}

static inline void evindex_s_set(struct rbindex *evidx, struct event *e)
{
	rbindex_s_set(evidx, e);
}

void evindex_s_seek(struct rbindex *evidx, double t);

static inline void evindex_s_rewind(struct rbindex *evidx)
{
	evidx->cur_s = evidx->ls.front->next;
//	evidx->cur_s = evidx->ls.front;
}

static inline int evindex_s_final(struct rbindex *evidx)
{
	evidx->cur_s = evidx->rsentinel;
}

static inline int evindex_s_end(struct rbindex *evidx)
{
	return evidx->cur_s == NULL;
}

static inline void evindex_seq_on(struct rbindex *evidx)
{
	rbindex_seq_on(evidx);
}

void evindex_seq_off(struct rbindex *evidx);

static inline struct event *evindex_find(rb_traverser *it, struct rbindex *evidx, struct event *key)
{
	return rbindex_find(it, evidx, key);
}

static inline void evindex_s_delete(struct rbindex *evidx, struct event *e)
{
	rbindex_s_delete(evidx, (void *)e);
}

static inline void evindex_rb_delete(struct rbindex *evidx, struct event *e)
{
	rbindex_rb_delete(evidx, (void *)e);
}

static inline void evindex_delete(struct rbindex *evidx, struct event *e)
{
	rbindex_delete(evidx, (void *)e);
}

static inline void evindex_rb_insert(struct rbindex *evidx, struct event *e)
{
	rbindex_rb_insert(evidx, (void *)e);
}

static inline void evindex_s_insert(struct rbindex *evidx, struct event *e)
{
	rbindex_s_insert(evidx, (void *)e);
}

static inline void evindex_insert(struct rbindex *evidx, struct event *e)
{
	if(rbindex_isseq(evidx))
		evindex_s_insert(evidx, e);
	else
		evindex_rb_insert(evidx, e);
}

void evindex_destroy(struct genealogy *G, struct rbindex *evidx);
struct rbindex *evindex_create(struct genealogy *G, struct config *cfg);

#endif
