#ifndef EINDEX_H
#define EINDEX_H

#include "rbindex.h"
#include "smc.h"
#include "global.h"

static inline struct edge *eindex_cur(seq_traverser it)
{
	return (struct edge *)rbindex_cur(it);
}

static inline struct edge *eindex_prev(seq_traverser *it)
{
	return (struct edge *)rbindex_prev(it);
}

static inline struct edge *eindex_next(seq_traverser *it)
{
	return (struct edge *)rbindex_next(it);
}

static inline void eindex_s_set(struct rbindex *eidx, struct edge *e)
{
	rbindex_s_set(eidx, e);
}
void eindex_s_seek_ttop(struct rbindex *eidx, double ttop);
void eindex_s_seek_t(struct rbindex *eidx, double ttop, double tbot);
void eindex_s_seek(struct rbindex *eidx, double ttop, double tbot, int eid);

static inline void eindex_s_rewind(struct rbindex *eidx)
{
	eidx->cur_s = eidx->ls.front->next;
}

static inline int eindex_s_final(struct rbindex *eidx)
{
	eidx->cur_s = eidx->rsentinel;
}

static inline void eindex_seq_on(struct rbindex *eidx)
{
	rbindex_seq_on(eidx);
}

void eindex_seq_off(struct rbindex *eidx);

static inline struct edge *eindex_find(seq_traverser *it, struct rbindex *eidx, struct edge *key)
{
	return rbindex_find(it, eidx, key);
}

static inline void eindex_s_delete(struct rbindex *eidx, struct edge *e)
{
	rbindex_s_delete(eidx, (void *)e);
}

static inline void eindex_rb_delete(struct rbindex *eidx, struct edge *e)
{
	rbindex_rb_delete(eidx, (void *)e);
}

static inline void eindex_delete(struct rbindex *eidx, struct edge *e)
{
	rbindex_delete(eidx, (void *)e);
}

static inline void eindex_rb_insert(struct rbindex *eidx, struct edge *e)
{
	rbindex_rb_insert(eidx, (void *)e);
}

void eindex_s_insert(struct rbindex *eidx, struct edge *e);

static inline void eindex_insert(struct rbindex *eidx, struct edge *e)
{
	if(rbindex_isseq(eidx))
		eindex_s_insert(eidx, e);
	else
		eindex_rb_insert(eidx, e);
}

void eindex_reset(struct genealogy *G, struct rbindex *eidx);
void eindex_destroy(struct genealogy *G, struct rbindex *eidx);
struct rbindex *eindex_create(struct genealogy *G, int pop);

#endif
