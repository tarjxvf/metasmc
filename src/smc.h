#ifndef SIMULATOR_H
#define SIMULATOR_H

#include "global.h"
#include "mutation.h"
#include "rbindex.h"
#include "tsindex.h"
#include "evindex.h"

#define NODE_COAL	0
#define NODE_MIGR	1
#define NODE_XOVER	2
#define NODE_SAM	3
#define NODE_FLOAT	4
#define NODE_DUMMY	5

#define AS_COAL_NODE(n)	((struct coal_node *)(n))
#define AS_MIGR_NODE(n)	((struct migr_node *)(n))
#define AS_XOVER_NODE(n)	((struct xover_node *)(n))
#define AS_SAM_NODE(n)	((struct sam_node *)(n))
#define AS_FLOAT_NODE(n)	((struct dummy_node *)(n))
#define AS_DUMMY_NODE(n)	((struct dummy_node *)(n))

static inline int iscoalnode(struct node *nd)
{
	return nd->type == NODE_COAL;
}

static inline int isxovernode(struct node *nd)
{
	return nd->type == NODE_XOVER;
}

static inline int ismigrnode(struct node *nd)
{
	return nd->type == NODE_MIGR;
}

static inline int issamnode(struct node *nd)
{
	return nd->type == NODE_SAM;
}

static inline int isfloatnode(struct node *nd)
{
	return nd->type == NODE_FLOAT;
}

static inline int isdummynode(struct node *nd)
{
	return nd->type == NODE_DUMMY;
}

static int node_flag_gettype(struct node *nd)
{
	return nd->type;
}

static void edge_flag_settype(struct node *nd, char type)
{
	nd->type = type;
}

static inline char isdeleted(struct node *e)
{
	return e->deleted;
}

static inline void edge_flag_delete(struct node *e)
{
	e->deleted = 1;
}

static inline void edge_flag_undelete(struct node *e)
{
	e->deleted = 0;
}

static inline void edge_flag_setdeleted(struct node *e, char flag)
{
	e->deleted = flag;
}

static inline int edge_flag_getitop(struct node *e)
{
	return e->itop;
}

static inline void edge_flag_setleft(struct node *e)
{
	e->itop = 0;
}

static inline void edge_flag_setright(struct node *e)
{
	e->itop = 1;
}

static inline void edge_flag_setitop(struct node *e, char flag)
{
	e->itop = flag;
}

void free_node(struct genealogy *G, struct node *nd);
struct node *alloc_node(struct genealogy *G, int type, int pop, double t);

struct node_set {
	int maxn;
	int n;
	struct node **nodes;
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
	struct node **eptrs;		// Array of node pointers
	struct list idx_queue;	// Queue of index in eptrs
	struct list id_list;
	struct list e_delete_list;
};

struct genealogy {
	int nsam;
	struct config *cfg;
	struct population *pops;
//	struct list *evlist;	// List of events
	struct evindex *evidx;
	struct event *ev_dxvr;	// Dummy recombination event occuring above localMRCA and below root.
	int edgeid;

	int nnodes;		// Number of nodes in local tree
	struct list n_delete_list;

//	struct list_head *evlcurr;
	struct node *root;
	double troot;		// height of existing tree
	double t;
	double total;		// Total length of the local tree

//	struct node **pTreeEdgesToCoalesceArray;
	struct node *localMRCA;
	int lb;
	int ub;
	int x;

	struct tsindex *tr_xover;
	struct node_set *trunk;

	int maxR;
	int nR[2];
	int *R[2];
	int curridx;
};

static inline void insert_event_join_increase(struct genealogy *G, struct join_event *jev)
{
	jev->dn[jev->popi]++;
	jev->dn[jev->popj]--;
}

static inline void insert_event_splt_increase(struct genealogy *G, struct splt_event *sev)
{
	sev->dn[sev->pop]--;
	sev->dn[sev->newpop]++;
}

void insert_event_rb_join(struct genealogy *G, struct join_event *jev);
void insert_event_rb_splt(struct genealogy *G, struct splt_event *sev);

static inline void insert_event_join(struct genealogy *G, struct event *ev)
{
	if(!rbindex_isseq(G->evidx->idx))
		insert_event_rb_join(G, (struct join_event *)ev);
	else
		insert_event_join_increase(G, (struct join_event *)ev);
}

static inline void insert_event_splt(struct genealogy *G, struct event *ev)
{
	if(!rbindex_isseq(G->evidx->idx))
		insert_event_rb_splt(G, (struct splt_event *)ev);
	else
		insert_event_splt_increase(G, (struct splt_event *)ev);
}

static inline void insert_event_rb(struct genealogy *G, struct event *ev)
{
//	if(ev->type == EVENT_JOIN){
//		insert_event_rb_join(G, ev);

//	}else if(ev->type == EVENT_SPLT){
//		insert_event_rb_splt(G, ev);

//	}else{
		evindex_rb_insert(G->evidx, ev);
//	}
}

static inline void insert_event(struct genealogy *G, struct event *ev)
{
//	if(ev->type == EVENT_JOIN){
//		insert_event_join(G, ev);

//	}else if(ev->type == EVENT_SPLT){
//		insert_event_splt(G, ev);

//	}else{
		evindex_insert(G->evidx, ev);
//	}
}

static inline void remove_event_join_decrease(struct genealogy *G, struct join_event *jev)
{
	jev->dn[jev->popi]--;
	jev->dn[jev->popj]++;
}

static inline void remove_event_splt_decrease(struct genealogy *G, struct splt_event *sev)
{
	sev->dn[sev->pop]++;
	sev->dn[sev->newpop]--;
}

void remove_event_rb_join(struct genealogy *G, struct join_event *jev);
void remove_event_rb_splt(struct genealogy *G, struct splt_event *sev);

static inline void remove_event_rb(struct genealogy *G, struct event *ev)
{
	evindex_rb_delete(G->evidx, ev);
	free_event(G->cfg, ev);
//	free(GET_LIST(ev));
}

static inline void remove_event_s(struct genealogy *G, struct event *ev)
{
	evindex_s_delete(G->evidx, ev);
	free_event(G->cfg, ev);
//	free(GET_LIST(ev));
}

static inline void remove_event(struct genealogy *G, struct event *ev)
{
	evindex_delete(G->evidx, ev);
	free_event(G->cfg, ev);
//	free(GET_LIST(ev));
}

static inline void remove_event_rb_josp(struct genealogy *G, struct event *ev)
{
	if(ev->type == EVENT_JOIN)
		remove_event_rb_join(G, (struct join_event *)ev);
	else if(ev->type == EVENT_SPLT)
		remove_event_rb_splt(G, (struct splt_event *)ev);
	else
		remove_event_rb(G, ev);
}

static inline void remove_event_s_josp(struct genealogy *G, struct event *ev)
{
	if(ev->type == EVENT_JOIN)
		remove_event_join_decrease(G, (struct join_event *)ev);
	else if(ev->type == EVENT_SPLT)
		remove_event_splt_decrease(G, (struct splt_event *)ev);
	else
		remove_event_s(G, ev);
}

static inline void remove_event_josp(struct genealogy *G, struct event *ev)
{
	if(!rbindex_isseq(G->evidx->idx))
		remove_event_rb_josp(G, ev);
	else
		remove_event_s_josp(G, ev);
}

//int simulate(struct reference *, struct genealogy *, int, struct frag *);
int simulate(struct genealogy *G, struct profile *prof);

struct genealogy *alloc_genealogy(struct config *, struct profile *);
void destroy_genealogy(struct genealogy *);
void clear_genealogy(struct genealogy *);

struct event *rnd_select_point(struct genealogy *G, struct node **eo, int *popo, double *to);

#endif
