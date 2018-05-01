#ifndef MUTATION_H
#define MUTATION_H

#include "global.h"
#include "smc.h"

#define MODEL_JC69	1
#define MODEL_K80	2
#define MODEL_F81	3
#define MODEL_HKY85	4
#define MODEL_T92	5
#define MODEL_TN93	6
#define MODEL_GTR	7

struct mutation {
	/*** Mutation model ***/
	int model;					// Mutation model
	double theta;
	double mpar[SQNUCS - NUM_NUCS];			// Parameters of mutation model
	double pi[NUM_NUCS];				// Equalibrium allele frequencies. All entries are 0.25 by default
	double D[NUM_NUCS];				// Eigenvalues of mutation matrix
	double Cijk[NUM_NUCS][NUM_NUCS][NUM_NUCS];	// Coefficient term in matrix exponential
};

extern void (*generate_sequence)(struct reference *, struct genealogy *, int, int);
void generate_sequence_infinite_fast(struct reference *, struct genealogy *, int, int);
void generate_sequence_infinite_slow(struct reference *, struct genealogy *, int, int);

void generate_sequence_model_fast(struct reference *, struct genealogy *, int, int);
void generate_sequence_model_slow(struct reference *, struct genealogy *, int, int);

void init_mutation_model(struct mutation *mmut);
void dump_mutation_model(struct mutation *mmut);
void setmatrix(struct mutation *mmut, double *Pmut, double t);

extern char *model_name[];// Model names
extern int npar[];	// Length of parameter list of all models

#endif
