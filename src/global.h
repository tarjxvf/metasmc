#ifndef GLOBAL_H
#define GLOBAL_H

#include <stdio.h>
#include "slab.h"

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

//#define EVENT_MMUT	10	/* Change of mutation model */

struct config;
struct genealogy;
struct mutation;

struct event {
	int type;
	double t;
};

struct coal_event {
	int type;	// type==EVENT_COAL
	double t;
	int pop;
};

struct migr_event {
	int type;	// type==EVENT_MIGR
	double t;
	int spop;
	int dpop;
};

struct grow_event {
	int type;	// type==EVENT_GROW
	double t;
	int pop;
	double alpha;
};

struct size_event {
	int type;	// type==EVENT_SIZE
	double t;
	int pop;
	double size;
};

struct gmig_event {
	int type;
	double t;
	double rmig;
};

struct rmig_event {
	int type;
	double t;
	int popi;
	int popj;
	double rmig;
};

struct gsiz_event {
	int type;
	double t;
	double size;
};

struct ggro_event {
	int type;
	double t;
	double alpha;	// Growth rate of all subpopulations
};

struct join_event {
	int type;
	double t;
	int popi;	// Subpopulation to be absorbed
	int popj;
};

struct splt_event {
	int type;
	double t;
	int pop;	// Subpopulation to be splitted
	int newpop;	// New subpopulation
	double prop;	// proportion is probability that each lineage stays in pop-i. (p, 1-p are admixt. proport.
};

/*struct mmut_event {
	int type;
	double t;
	int pop;
	struct mutation *mmut;
};*/
void remove_event(struct genealogy *, struct event *);
//struct event *alloc_event(struct genealogy *, int, double);
struct event *alloc_event(struct config *, int, double);
void print_event(struct event *ev);

struct list_head {
	struct list_head *next;
	struct list_head **prev;
};

struct list {
	struct list_head *front;
	struct list_head **rear;
	int n;
};

void list_init(struct list *ls);
void list_concat(struct list *dst, struct list *src);

/* List operations using object pointer. */
void list_insbefore(struct list_head *ref, void *item);
void list_remove(struct list *ls, void *item);
void list_add(struct list *ls, void *item);	// Add an item at list head
void list_append(struct list *ls, void *item);	// Add an item at list end

/* List operations using list_head pointer. */
static __inline__ void __list_add(struct list *ls, struct list_head *l) 
{
	if(ls->front)
		ls->front->prev = &l->next;
	else
		ls->rear = &l->next;

	l->next = ls->front;
	ls->front = l;
	l->prev = &ls->front;
	ls->n++;
}

/* Insert an item before an item. */
static __inline__ void __list_insbefore(struct list_head *ref, struct list_head *l)
{
	l->next = ref;
	l->prev = ref->prev;
	ref->prev = &l->next;
	*l->prev = l;
}

/* Append an item after a  list */
static __inline__ void __list_append(struct list *ls, struct list_head *l)
{
	l->next = NULL;
	*ls->rear = l;
	l->prev = ls->rear;
	ls->rear = &l->next;
	ls->n++;
}

/* Remove an item from the list. */
static __inline__ void __list_remove(struct list *ls, struct list_head *l)
{
	*l->prev = l->next;
	if(l->next)
		l->next->prev = l->prev;
	else	// If the removed item is the last one
		ls->rear = l->prev;
	ls->n--;
}

/*void list_insbefore(struct list_head *ref, struct list_head *item);
void list_remove(struct list_head *item);
void list_add(struct list_head **head, struct list_head *item);
void list_append(struct list_head **head, struct list_head *item);*/

#define GET_LIST(obj)	((struct list_head *)((char *)(obj) - sizeof(struct list_head)))
#define GET_OBJ(l)	((char *)(l) + sizeof(struct list_head))

/* Single-end read. */
struct read{
	int start;
	int end;
	char *seq;
	char *qual;
};

/* Fragment of paired-end reads. */
struct frag{
	int id;
//	int chr;	// Chromosome number
//	double start;
//	double end;
	int start;	// Chromosome start position
	int end;
	int pop;
	int nread;
	struct read *rd;
	int trunk;
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

#define NUM_NUCS	4
#define SQNUCS	NUM_NUCS * NUM_NUCS

struct profile {
	int npop;
//	double *popsize;

	char *reffile;
//	int nchr;
//	int reflen;
	struct reference *ref;
	int chrnum;	// Index of chromosome in the reference;

	int *nfrag_pop;
	int nfrag;
	struct frag *fgset;
};

int fgcompar(const void *a, const void *b);
void print_profile(struct profile *prof, FILE *outfp);
struct profile *generate_profile(char *reffile, int chrnum, int npop, int *nfrags, int fraglen, int paired, int rdlen);
struct profile *load_profile(FILE *filp);
void unload_profile(struct profile *prof);

struct mutation;
struct config {
	unsigned int seed;
	int print_tree;	// 1 if you wish to print trees
	int gensam;	// 1 if you wish to generate sequences
	FILE *treefp;	// If print_tree is set to 1, this points to file object of tree output
	FILE *readfp;	// Output file for simulated reads

	/*** Basic parameters. ***/
	double rho;
	double tdiv;	// Divergence time

	/*** Demographic model. ***/
	int npop;	// Number of initial subpopulations
//	int *nsams;	// Number of samples in each population
	int nsplt;	// Number of subpopulations created by splitting events.
	int npop_all;	// Maximum number of populations, including those created by splitting events

	double *size;	// Initial subpopulation size (at time 0)
	double **mmig;	// Initial migration matrix(at time 0)
	double *grate;	// Initial growth rates

	/*** Mutation model ***/
	struct mutation **mmut;

	struct list evlist;

	/*** Caches of frequently-used objects (slab allocators) ***/
	struct mem_cache *node_cache[5];
	struct mem_cache *event_cache[2];
	struct mem_cache *edge_cache;
	struct mem_cache *frag_cache;

	int maxfrag;
int debug;
};

struct config *create_config(int seed, int print_tree, int gensam, FILE *, FILE *, int maxfrag, struct profile *prof, double rho);
void destroy_config(struct config *cfg);
void dump_config(struct config *cfg);

int register_mutation_model(struct config *cfg, int pop, struct mutation *mmut);

void print_fragment(FILE *outfp, struct frag *fg);	

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
//int add_event_mmut(struct config *cfg, double t, int pop, double theta, int model, double *pars, double *pi);

extern char nucl[];
int nucl_index(int ch);

double dunif01();
double dexp();
unsigned int dunif(unsigned int max);
int poisso(double u);
void seedit( char *flag );
//void init_rand();
void init_rand(unsigned int s);
void finish_rand();
void seed();

#endif
