#ifndef GLOBAL_H
#define GLOBAL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "list.h"
#include "cache.h"

#define NUM_THREADS 4

#define EVENT_COAL	0	/* Coalescent */
#define EVENT_MIGR	1	/* Migration */
#define EVENT_GROW	2	/* Change growth rate of a subpopulation. */
#define EVENT_SIZE	3	/* Change size of a subpopulation */
#define EVENT_RMIG	4	/* Change migration rate. */
#define EVENT_GMIG	5	/* Change migration rate. */
#define EVENT_GSIZ	6	/* Change size of all subpopulations. */
#define EVENT_GGRO	7	/* Change growth rate of all subpopulations. */
#define EVENT_JOIN	8	/* Population join */
#define EVENT_SPLT	9	/* Population split */
#define EVENT_DUMY	10	/* Old MRCA. */
#define EVENT_DXVR	11	/* Recombination on dummy lineages. */
#define EVENT_SAMP	12	/* Add new sample. */

struct config;
struct genealogy;
struct mutation;

struct coal_node;
struct node;

struct node_set {
	int maxn;
	int n;
	struct node **nodes;
};

struct event {
	// dn: Change of the number of lineages in affected populations. sumdn: Sum of dn values of the (red-black tree) nodes below this event
	double t;
	unsigned char type;
	unsigned char dn_off;
	unsigned char sumdn_off;
};

struct coal_event {
	// type==EVENT_COAL,  dn[pop] == -1 and dn[i] == 0 for other i
	double t;
	unsigned char type;
	unsigned char dn_off;
	unsigned char sumdn_off;

	int pop;
	struct coal_node *nd;
};

struct migr_event {
	// type==EVENT_MIGR, dn[dpop] == +1 and dn[spop] == -1
	double t;
	unsigned char type;
	unsigned char dn_off;
	unsigned char sumdn_off;

	int spop;
	int dpop;
	struct migr_node *nd;
};

struct grow_event {
	// type==EVENT_GROW and dn Must be zero
	double t;
	unsigned char type;
	unsigned char dn_off;
	unsigned char sumdn_off;

	int pop;
	double alpha;
};

struct size_event {
	// type==EVENT_SIZE and dn must be zero
	double t;
	unsigned char type;
	unsigned char dn_off;
	unsigned char sumdn_off;

	int pop;
	double size;
};

struct gmig_event {
	// dn must be zero
	double t;
	unsigned char type;
	unsigned char dn_off;
	unsigned char sumdn_off;

	double rmig;
};

struct rmig_event {
	// dn must be zero
	double t;
	unsigned char type;
	unsigned char dn_off;
	unsigned char sumdn_off;

	int popi;
	int popj;
	double rmig;
};

struct gsiz_event {
	// dn must be zero
	double t;
	unsigned char type;
	unsigned char dn_off;
	unsigned char sumdn_off;

	double size;
};

struct ggro_event {
	// dn must be zero
	double t;
	unsigned char type;
	unsigned char dn_off;
	unsigned char sumdn_off;

	double alpha;	// Growth rate of all subpopulations
};

struct join_event {
	// dn[0] depends on the number of lineages in source population. dn[1] = -dn[0]
	double t;
	unsigned char type;
	unsigned char dn_off;
	unsigned char sumdn_off;

	int popi;	// Subpopulation to be absorbed
	int popj;
	struct node_set ndls;
};

struct splt_event {
	// dn[0] depend on the number of migrated lineages. dn[1] = -dn[0]
	double t;
	unsigned char type;
	unsigned char dn_off;
	unsigned char sumdn_off;

	int pop;	// Subpopulation to be splitted
	int newpop;	// New subpopulation
	double prop;	// proportion is probability that each lineage stays in pop-i. (p, 1-p are admixt. proport.
	struct node_set ndls;
};

struct samp_event {
	// Equals to the number of samples of each population added at time t
	double t;
	unsigned char type;
	unsigned char dn_off;
	unsigned char sumdn_off;
};

#define GET_DN(ev) ((int *)((char *)(ev) + ((struct event *)(ev))->dn_off))
#define GET_SUMDN(ev) ((int *)((char *)(ev) + ((struct event *)(ev))->sumdn_off))

struct event *alloc_event(struct config *, int, double);
void free_event(struct config *cfg, struct event *ev);
void print_event(struct config *cfg, struct event *ev);

/* Single-end read. */
struct read{
	int start;
	int end;
	char *seq;
	char *qual;
};

struct fginfo{
	int end;
	unsigned int pop:8;
	unsigned int nread:8;
	unsigned int trunk:1;
};

struct reference {
	FILE *filp;
	int nchr;
	int *seq_start;	// Start position of each chromosome
	int *chrlen;
	int currchr;
	int curr;
};

struct reference *load_reference(char *file);
void unload_reference(struct reference *ref);
void reload_reference(struct reference *ref, int chrnum);
void load_chr(struct reference *ref, int chrnum, char **strp);

#define NUM_NUCS	4
#define SQNUCS	NUM_NUCS * NUM_NUCS

struct profile {
	int npop;

	char *reffile;
	struct reference *ref;
	int chrnum;	// Index of chromosome in the reference;

	int *nfrag_pop;
	int *ntrunks;
	int nfrag;

	int *fgstart;
//	int *fgend;
	int *fgid;
	struct fginfo *info;
	struct sam_node **nds;
	struct read **rdset;
};

int fgcompar(const void *a, const void *b);
void print_profile(struct profile *prof, FILE *outfp);
struct profile *generate_profile(char *reffile, int chrnum, int npop, int *nfrags, int fraglen, int paired, int rdlen, int *ntrunks);
struct profile *load_profile(FILE *filp, int genseq);
void unload_profile(struct profile *prof);

struct mutation;
struct config {
	struct profile *prof;
	unsigned char seed;
	int print_tree;	// 1 if you wish to print trees
	int gensam;	// 1 if you wish to generate sequences
	FILE *treefp;	// If print_tree is set to 1, this points to file object of tree output
	FILE *readfp;	// Output file for simulated reads

	/*** Basic parameters. ***/
	double rho;
	double tdiv;	// Divergence time

	/*** Demographic model. ***/
	int npop;	// Number of initial subpopulations
	int nsplt;	// Number of subpopulations created by splitting events.
	int npop_all;	// Maximum number of populations, including those created by splitting events

	double *size;	// Initial subpopulation size (at time 0)
	double **mmig;	// Initial migration matrix(at time 0)
	double *grate;	// Initial growth rates

	/*** Mutation model ***/
	struct mutation **mmut;

	struct list evlist;
	int ndevents;		// Number of demographic events
	struct event **devents;	// Array of demographic events

	/*** Caches of frequently-used objects ***/
	struct cache *node_cache[6];
	struct cache *event_cache[13];

	int maxfrag;
int debug;
};

struct config *create_config(int seed, int print_tree, int gensam, FILE *, FILE *, int maxfrag, struct profile *prof, double rho);
void destroy_config(struct config *cfg);
void dump_config(struct config *cfg);

int register_mutation_model(struct config *cfg, int pop, struct mutation *mmut);

//void print_fragment(FILE *outfp, struct frag *fg, struct read *rd);
void print_fragment(FILE *outfp, int fgstart, int fgid, struct fginfo *fgi, struct read *rdset);

int set_growth_rates(struct config *cfg, double *grate);
int set_migration_matrix(struct config *cfg, double *mmig);

int add_event_ggro(struct config *cfg, double t, double alpha);
int add_event_grow(struct config *cfg, double t, int pop, double alpha);
int add_event_gmig(struct config *cfg, double t, double rmig);
int add_event_rmig(struct config *cfg, double t, int popi, int popj, double rmig);
int add_event_gsiz(struct config *cfg, double t, double size);
int add_event_join(struct config *cfg, double t, int popi, int popj);
int add_event_splt(struct config *cfg, double t, int pop, double prop);
int add_event_size(struct config *cfg, double t, int pop, double size);
int add_event_samp(struct config *cfg, double t, int pop, double size);

extern char nucl[];
int nucl_index(int ch);

typedef size_t map_t;
#define CELLSIZE 32
#define NCELL_PER_MAP sizeof(map_t)

struct node {
	// type: NODE_COAL, NODE_MIGR, NODE_SAM, NODE_FLOAT
	double t;
	int set_id;
	int xtid;
	int idx;
	struct{
		unsigned char type:4;
		unsigned char itop:1;
		unsigned char deleted:1;
		unsigned char visited:2;
		unsigned char pop;
	};
	struct event *ev;
	struct node *in;
	struct node *out[2];
};

// node representing coalescent event
struct coal_node {
	// type==NODE_COAL
	double t;
	int set_id;
	int xtid;
	int idx;
	struct{
		unsigned char type:4;
		unsigned char itop:1;
		unsigned char deleted:1;
		unsigned char visited:2;
		unsigned char pop;
	};
	struct coal_event *ev;
	struct node *in;
	struct node *out[2];	//Edges below the node
	char *seq;
	map_t *mapped;
};

// Node representing recombination event. Not used in current implementation.
struct xover_node {
	// type==NODE_XOVER
	double t;
	int set_id;
	int xtid;
	int idx;
	struct{
		unsigned char type:4;
		unsigned char itop:1;
		unsigned char deleted:1;
		unsigned char visited:2;
		unsigned char pop;
	};
	struct event *ev;
	struct node *in_new;
	struct node *out;
	struct node *in;
};

// Note representing migration event
struct migr_node {
	// type==NODE_MIGR
	double t;
	int set_id;
	int xtid;
	int idx;
	struct{
		unsigned char type:4;
		unsigned char itop:1;
		unsigned char deleted:1;
		unsigned char visited:2;
		unsigned char pop;
	};
	struct migr_event *ev;
	struct node *in;
	struct node *out;
	int mgid;
};

// Node representing sample.
struct sam_node {
	// type==NODE_SAM
	double t;
	int set_id;
	int xtid;
	int idx;
	struct{
		unsigned char type:4;
		unsigned char itop:1;
		unsigned char deleted:1;
		unsigned char visited:2;
		unsigned char pop;
	};
	struct event *ev;
	struct node *in;
	int fgid;
};

// Node representing tip of dummy lineage which represents trapped ancestral material. Recombination is allowed on this type of lineage but take no effect.
struct dummy_node {
	// type==NODE_DUMMY of type==NODE_FLOAT
	double t;
	int set_id;
	int xtid;
	int idx;
	struct{
		unsigned char type:4;
		unsigned char itop:1;
		unsigned char deleted:1;
		unsigned char visited:2;
		unsigned char pop;
	};
	struct event *ev;
	struct node *in;
	struct node *out;
};

#define NODE_FLAG_VISITED_LEFT 0x1
#define NODE_FLAG_VISITED_RIGHT 0x2

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

static inline void migr_set_add(struct node_set *set, struct node *e)
{
	if(set->n >= set->maxn){
		set->maxn *= 2;
		set->nodes = realloc(set->nodes, sizeof(struct node *) * set->maxn);
	}
	AS_MIGR_NODE(e)->mgid = set->n;
	set->nodes[set->n++] = e;
}

static inline void migr_set_remove(struct node_set *set, struct node *e)
{
	int i;
	i = AS_MIGR_NODE(e)->mgid;
	set->nodes[i] = set->nodes[--(set->n)];
	AS_MIGR_NODE(set->nodes[i])->mgid = i;
}

static inline void node_set_init(struct node_set *set, int maxn)
{
	set->n = 0;
	set->maxn = maxn;
	set->nodes = malloc(sizeof(struct node *) * maxn);
}

static inline void node_set_resize(struct node_set *set, int maxn)
{
	if(maxn > set->maxn){
		set->maxn = maxn;
		set->nodes = realloc(set->nodes, sizeof(struct node *) * maxn);
//		memset(set->nodes, 0, sizeof(struct node *) * maxn);
	}
}

static inline void node_set_clear(struct node_set *set)
{
	set->n = 0;
}

static inline void node_set_destroy(struct node_set *set)
{
	set->maxn = set->n = 0;
	free(set->nodes);
}

static inline void node_set_add(struct node_set *set, struct node *e)
{
	e->set_id = set->n;
	set->nodes[set->n++] = e;
}

static inline struct node *node_set_get(struct node_set *set, int i)
{
	return set->nodes[i];
}

static inline void node_set_replace(struct node_set *set, int i, struct node *e)
{
	set->nodes[i]->set_id = -1;
	e->set_id = i;
	set->nodes[i] = e;
}

static inline struct node *node_set_remove(struct node_set *set, int i)
{
	struct node *e;

	e = set->nodes[i];
	set->nodes[i] = set->nodes[--(set->n)];
	set->nodes[i]->set_id = i;
//	set->nodes[set->n] = NULL;

	return e;
}

#endif
