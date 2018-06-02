#ifndef SIMULATOR_H
#define SIMULATOR_H

#include "global.h"
#include "mutation.h"
#include "rbindex.h"
#include "tsindex.h"

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

void free_node(struct genealogy *G, struct node *nd);
struct node *alloc_node(struct genealogy *G, int type, int pop, double t);
void free_edge(struct genealogy *G, struct edge *e);
struct edge *alloc_edge(struct genealogy *G, struct node *top, struct node *bot);

struct edge {
	struct node *top;
	struct node *bot;
	int itop;
	int eid;
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

struct genealogy {
	int nsam;
	struct config *cfg;
	struct population *pops;
	struct list *evlist;	// List of events
	struct event *ev_dxvr;	// Dummy recombination event occuring above localMRCA and below root.
	int edgeid;

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
