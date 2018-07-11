#ifndef GLOBAL_H
#define GLOBAL_H

#include <stdio.h>
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

//#define EVENT_MMUT	10	/* Change of mutation model */

//#define MAX(a, b) ((a) > (b))?(a):(b)

struct config;
struct genealogy;
struct mutation;

struct event {
	int type;
	double t;
	int *dn;	// Change of the number of lineages in affected populations
	int *sumdn;	// Sum of dn values of the (red-black tree) nodes below this event
};

struct coal_event {
	int type;	// type==EVENT_COAL
	double t;
	int *dn;	// dn[pop] == -1 and dn[i] == 0 for other i
	int *sumdn;
	int pop;
	struct coal_node *nd;
};

struct migr_event {
	int type;	// type==EVENT_MIGR
	double t;
	int *dn;	// dn[dpop] == +1 and dn[spop] == -1
	int *sumdn;
	int spop;
	int dpop;
	struct migr_node *nd;
};

struct grow_event {
	int type;	// type==EVENT_GROW
	double t;
	int *dn;	// Must be zero
	int *sumdn;
	int pop;
	double alpha;
};

struct size_event {
	int type;	// type==EVENT_SIZE
	double t;
	int *dn;	// Must be zero
	int *sumdn;
	int pop;
	double size;
};

struct gmig_event {
	int type;
	double t;
	int *dn;	// Must be zero
	int *sumdn;
	double rmig;
};

struct rmig_event {
	int type;
	double t;
	int *dn;	// Must be zero
	int *sumdn;
	int popi;
	int popj;
	double rmig;
};

struct gsiz_event {
	int type;
	double t;
	int *dn;	// Must be zero
	int *sumdn;
	double size;
};

struct ggro_event {
	int type;
	double t;
	int *dn;	// Must be zero
	int *sumdn;
	double alpha;	// Growth rate of all subpopulations
};

struct join_event {
	int type;
	double t;
	int *dn;	// Depends on the number of lineages in source population. dn[1] = -dn[0]
	int *sumdn;
	int popi;	// Subpopulation to be absorbed
	int popj;
	struct list ndls;
};

struct splt_event {
	int type;
	double t;
	int *dn;	// Depend on the number of migrated lineages. dn[1] = -dn[0]
	int *sumdn;
	int pop;	// Subpopulation to be splitted
	int newpop;	// New subpopulation
	double prop;	// proportion is probability that each lineage stays in pop-i. (p, 1-p are admixt. proport.
	struct list ndls;
};

struct samp_event {
	int type;
	double t;
	int *dn;	// Equals to the number of samples of each population added at time t
	int *sumdn;
};

/*struct mmut_event {
	int type;
	double t;
	int pop;
	struct mutation *mmut;
};*/

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

struct sam_node;

/* Fragment of paired-end reads. */
struct frag{
	int start;	// Chromosome start position
	int end;
	unsigned int pop:8;
	unsigned int nread:8;
	unsigned int trunk:1;
	int id;
	struct sam_node *nd;
//	struct read *rd;
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

	struct frag *fgset;
	struct read **rdset;
};

int fgcompar(const void *a, const void *b);
void print_profile(struct profile *prof, FILE *outfp);
struct profile *generate_profile(char *reffile, int chrnum, int npop, int *nfrags, int fraglen, int paired, int rdlen, int *ntrunks);
struct profile *load_profile(FILE *filp);
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

void print_fragment(FILE *outfp, struct frag *fg, struct read *rd);

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

struct node;

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
//	struct frag *fg;	// pointer to corresponding fragment
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

#endif
