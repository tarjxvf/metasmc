#ifndef SIMULATOR_H
#define SIMULATOR_H

#include "global.h"
#include "mutation.h"
#include "rb.h"

#define NODE_COAL	0
#define NODE_MIGR	1
#define NODE_XOVER	2
#define NODE_SAM	3
#define NODE_FLOAT	4
#define NODE_DUMMY	5

struct node;
struct edge;

typedef size_t map_t;
#define CELLSIZE 32
#define NCELL_PER_MAP sizeof(map_t)

struct node {
	int type;	// type: NODE_COAL, NODE_MIGR, NODE_SAM, NODE_FLOAT
	double t;
	int pop;
	struct event *ev;
	struct edge *in;	// Edge above the node
};

// node representing coalescent event
struct coal_node {
	int type;	// type==NODE_COAL
	double t;
	int pop;
	struct coal_event *ev;
	struct edge *in;
	struct edge *out[2];	//Edges below the node
	char *seq;
	map_t *mapped;
};

// Node representing recombination event. Not used in current implementation.
struct xover_node {
	int type;	// type==NODE_XOVER
	double t;
	int pop;
	struct event *ev;
	struct edge *in_new;
	struct edge *out;
	struct edge *in;
};

// Note representing migration event
struct migr_node {
	int type;	// type==NODE_MIGR
	double t;
	int pop;
	struct migr_event *ev;
	struct edge *in;
	struct edge *out;
};

// Node representing sample.
struct sam_node {
	int type;	// type==NODE_SAM
	double t;
	int pop;
	struct event *ev;
	struct edge *in;
	struct frag *fg;	// pointer to corresponding fragment
};

// Node representing tip of dummy lineage which represents trapped ancestral material. Recombination is allowed on this type of lineage but take no effect.
struct dummy_node {
	int type;	// type==NODE_DUMMY of type==NODE_FLOAT
	double t;
	int pop;
	struct event *ev;
	struct edge *in;
	struct edge *out;
};

#define AS_COAL_NODE(n)	((struct coal_node *)(n))
#define AS_MIGR_NODE(n)	((struct migr_node *)(n))
#define AS_XOVER_NODE(n)	((struct xover_node *)(n))
#define AS_SAM_NODE(n)	((struct sam_node *)(n))
#define AS_FLOAT_NODE(n)	((struct dummy_node *)(n))
#define AS_DUMMY_NODE(n)	((struct dummy_node *)(n))

struct edge {
	struct node *top;
	struct node *bot;
	int itop;
	int xtid;	// Index of edge in binary indexed tree
	int idx;	// Index of the edge in eptr array of population
};

struct edge_set {
	int maxn;
	int n;
	struct edge **edges;
};

struct population {
	int n;
	double grate;			// Growth rate at tlast
	double size;			// Relative size at tlast
	double tlast;
	double *mrate;
	int enabled;			// Whether this population is enabled at tlast

	int nsam;			// Number of sample nodes

	struct mutation *mmut;

	/***** MaCS-like data structures. *****/
	int nedges;
	int maxedges;
	struct edge **eptrs;		// Array of edge pointers
	struct list idx_queue;	// Queue of index in eptrs

	/***** Red-black tree index of edges. The tree is ordered by times of top nodes. *****/
//	struct rb_table *etree;
	struct rbindex *eidx;
};

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

static inline int tsindex_dirty(struct tsindex *tr)
{
	return tr->flags & TSINDEX_DIRTY;
}

static inline int tsindex_isrebuild(struct tsindex *tr)
{
	return tr->flags & TSINDEX_REBUILD;
}

void tsindex_rebuild(struct tsindex *tr);
struct tsindex *tsindex_alloc(int nedges);
void tsindex_reset(struct tsindex *tr);
void tsindex_add(struct tsindex *tr, struct edge *e);
void tsindex_update(struct tsindex *tr, struct edge *e, double diff);
void tsindex_clear(struct tsindex *tr, struct edge *e);
void tsindex_free(struct tsindex *tr);

#define RBINDEX_BATCH	0x1	/* If this flag is set, the index is in batch mode. The red-black tree will not be built until rbindex_build is called. */
//typedef struct rb_traverser rb_traverser;
typedef struct list_head * rb_traverser;

struct rbindex {
	int flags;
	struct rb_table *tree;
	struct list ls;

	/* For batch mode. */
	int maxobj;
	int nobj;
	void **objs;
};

static inline void rbindex_setflag(struct rbindex *eidx, int flag)
{
	eidx->flags |= flag;
}

static inline void rbindex_clearflag(struct rbindex *eidx, int flag)
{
	eidx->flags &= ~flag;
}

static inline int rbindex_isbatch(struct rbindex *eidx)
{
	return eidx->flags & RBINDEX_BATCH;
}

struct genealogy {
	int nsam;
	struct config *cfg;
	struct population *pops;
	struct list *evlist;	// List of events
	struct event *ev_dxvr;	// Dummy recombination event occuring above localMRCA and below root.

	int nedges;		// Number of edges in local tree
//	struct list e_list;	// List of edges in the local tree
	struct list n_list;	// List of sample nodes

	struct list_head *evlcurr;
	struct node *root;
	double troot;		// height of existing tree
	double t;
	double total;		// Total length of the local tree

	struct edge **pTreeEdgesToCoalesceArray;
	struct node *localMRCA;

	struct tsindex *tr_xover;
};

//int simulate(struct reference *, struct genealogy *, int, struct frag *);
int simulate(struct genealogy *G, struct profile *prof);

struct genealogy *alloc_genealogy(struct config *, struct profile *);
void destroy_genealogy(struct genealogy *);
void clear_genealogy(struct genealogy *);

struct event *rnd_select_point(struct genealogy *G, struct edge **eo, int *popo, double *to);

#endif
