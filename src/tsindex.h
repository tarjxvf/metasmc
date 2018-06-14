#ifndef TSINDEX_H
#define TSINDEX_H

#include "bit.h"
#include "global.h"

/* Tree size index. */
#define TSINDEX_REBUILD	0x1	// If this flag is set, the index is waiting for batch rebuild and any operations will not update binary index tree
#define TSINDEX_DIRTY	0x2	// This flag indicates that binary index tree is inconsistent to object list. The inconsistency must be fixed by tsindex_rebuild. Search operation is disabled if this flag is set.
struct tsindex {
	int flags;
	struct bit *index;
	int nedges;	// Equals to the number of occupied nodes in binary indexed tree
	int maxnodes;	// Equals to n in binary indexed tree
	int maxedges;
	struct edge **edges;
	struct list free_list;
};

static inline void tsindex_setflag(struct tsindex *tr, int flag)
{
	tr->flags |= flag;
}

static inline void tsindex_clearflag(struct tsindex *tr, int flag)
{
	tr->flags &= ~flag;
}

static inline int tsindex_isin(struct tsindex *tr, struct edge *e)
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

static inline double tsindex_size(struct tsindex *tr)
{
	if(tsindex_isrebuild(tr))
		return -1;
	else
		return bit_total(tr->index);
}

void tsindex_rebuild(struct tsindex *tr);
struct tsindex *tsindex_alloc(int nedges);
void tsindex_reset(struct tsindex *tr);
void tsindex_add(struct tsindex *tr, struct edge *e);
void tsindex_update(struct tsindex *tr, struct edge *e, double diff);
void tsindex_clear(struct tsindex *tr, struct edge *e);
void tsindex_free(struct tsindex *tr);
void tsindex_dump(struct tsindex *tr);
struct edge *tsindex_search(struct tsindex *tr, double g, double *cum);

#endif

