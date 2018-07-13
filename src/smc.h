#ifndef SIMULATOR_H
#define SIMULATOR_H

#include "global.h"
#include "mutation.h"
#include "rbindex.h"
#include "tsindex.h"
#include "evindex.h"

struct population {
	int n;
	double grate;			// Growth rate at tlast
	double size;			// Relative size at tlast
	double tlast;
	double *mrate;
	int enabled;			// Whether this population is enabled at tlast

	int nsam;			// Number of sample nodes

	struct mutation *mmut;

	int nedges;
	int maxedges;
	struct node **eptrs;		// Array of node pointers
	struct list idx_queue;	// Queue of index in eptrs
	struct list id_list;
};

struct genealogy {
	struct config *cfg;
	struct population *pops;
	struct evindex *evidx;
	struct event *ev_dxvr;	// Dummy recombination event occuring above localMRCA and below root.

	struct node *root;
	double troot;		// height of existing tree
	double t;
	double total;		// Total length of the local tree

	struct node *localMRCA;
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
	GET_DN(jev)[jev->popi]++;
	GET_DN(jev)[jev->popj]--;
}

static inline void insert_event_splt_increase(struct genealogy *G, struct splt_event *sev)
{
	GET_DN(sev)[sev->pop]--;
	GET_DN(sev)[sev->newpop]++;
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
	evindex_rb_insert(G->evidx, ev);
}

static inline void insert_event(struct genealogy *G, struct event *ev)
{
	evindex_insert(G->evidx, ev);
}

static inline void remove_event_join_decrease(struct genealogy *G, struct join_event *jev)
{
	GET_DN(jev)[jev->popi]--;
	GET_DN(jev)[jev->popj]++;
}

static inline void remove_event_splt_decrease(struct genealogy *G, struct splt_event *sev)
{
	GET_DN(sev)[sev->pop]++;
	GET_DN(sev)[sev->newpop]--;
}

void remove_event_rb_join(struct genealogy *G, struct join_event *jev);
void remove_event_rb_splt(struct genealogy *G, struct splt_event *sev);

static inline void remove_event_rb(struct genealogy *G, struct event *ev)
{
	evindex_rb_delete(G->evidx, ev);
	free_event(G->cfg, ev);
}

static inline void __remove_event_s(struct genealogy *G, struct event *ev)
{
	__evindex_s_delete(G->evidx, ev);
	free_event(G->cfg, ev);
}

static inline void remove_event_s(struct genealogy *G, struct event *ev)
{
	evindex_s_delete(G->evidx, ev);
	free_event(G->cfg, ev);
}

static inline void remove_event(struct genealogy *G, struct event *ev)
{
	evindex_delete(G->evidx, ev);
	free_event(G->cfg, ev);
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
