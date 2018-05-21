#ifndef SIMULATOR_H
#define SIMULATOR_H

#include "global.h"
#include "mutation.h"

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
};

/* Tree size index. */
struct tsindex {
	struct bit *index;
	int nedges;	// Equals to the number of occupied nodes in binary indexed tree
	int maxnodes;	// Equals to n in binary indexed tree
	int maxedges;
	struct edge **edges;
	struct list free_list;
};

struct tsindex *tsindex_alloc(int nedges);
void tsindex_reset(struct tsindex *tr);
void tsindex_add(struct tsindex *tr, struct edge *e);
void tsindex_update(struct tsindex *tr, struct edge *e, double diff);
void tsindex_clear(struct tsindex *tr, struct edge *e);
void tsindex_free(struct tsindex *tr);

struct genealogy {
	int nsam;
	struct config *cfg;
	struct population *pops;
	struct list *evlist;	// List of events
	struct event *ev_dxvr;	// Dummy recombination event occuring above localMRCA and below root.

	int nedges;		// Number of edges in local tree
	struct list e_list;	// List of edges in the local tree
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
