#ifndef TSINDEX_H
#define TSINDEX_H

#include "bit.h"
#include "global.h"

extern unsigned long long t_ts_rebuild;
extern unsigned long long n_ts_resize;

/* Tree size index. */
#define TSINDEX_REBUILD	0x1	// If this flag is set, the index is waiting for batch rebuild and any operations will not update binary index tree
#define TSINDEX_DIRTY	0x2	// This flag indicates that binary index tree is inconsistent to object list. The inconsistency must be fixed by tsindex_rebuild. Search operation is disabled if this flag is set.
struct tsindex {
	int flags;
	struct bit *index;
	int nedges;	// Equals to the number of occupied nodes in binary indexed tree
	int maxnodes;	// Equals to n in binary indexed tree
	int maxedges;
	struct node **edges;
	double *weights;
	struct list free_list;
	struct list id_list;
};

static inline void tsindex_setflag(struct tsindex *tr, int flag)
{
	tr->flags |= flag;
}

static inline void tsindex_clearflag(struct tsindex *tr, int flag)
{
	tr->flags &= ~flag;
}

static inline int tsindex_isin(struct tsindex *tr, struct node *e)
{
	return e->xtid > 0;
}

static inline int tsindex_dirty(struct tsindex *tr)
{
	return tr->flags & TSINDEX_DIRTY;
}

static inline int tsindex_isrebuild(struct tsindex *tr)
{
	return tr->flags & TSINDEX_REBUILD;
}

static inline struct node *tsindex_search(struct tsindex *tr, double g, double *cum)
{
	int eid;

	if(tsindex_isrebuild(tr)){
		return NULL;

	}else{
		eid = bit_getindex(tr->index, g);
		*cum = bit_cumfreq(tr->index, eid - 1);

		return tr->edges[eid];
	}
}

static inline void tsindex_update_r(struct tsindex *tr, struct node *e, double diff)
{
	tr->weights[e->xtid] += diff;
}

static inline void tsindex_update(struct tsindex *tr, struct node *e, double diff)
{
	if(!tsindex_isrebuild(tr))
		bit_update(tr->index, e->xtid, diff);
	tr->weights[e->xtid] += diff;
}

double tsindex_size(struct tsindex *tr);
void tsindex_rebuild(struct tsindex *tr);
struct tsindex *tsindex_alloc(int nedges);
void tsindex_reset(struct tsindex *tr);
void tsindex_replace_r(struct tsindex *tr, int id, struct node *e);
void tsindex_replace(struct tsindex *tr, int id, struct node *e);
void tsindex_add_r(struct tsindex *tr, struct node *e);
void tsindex_add(struct tsindex *tr, struct node *e);
void tsindex_clear_r(struct tsindex *tr, struct node *e);
void tsindex_clear(struct tsindex *tr, struct node *e);
void tsindex_free(struct tsindex *tr);
void tsindex_dump(struct tsindex *tr);
struct node *tsindex_search(struct tsindex *tr, double g, double *cum);

#endif

