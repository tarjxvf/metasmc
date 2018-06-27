#ifndef BIT_H
#define BIT_H

#include <stdio.h>

struct bit {
	int maxnode;
	int n;
	double *freq;
};

void __bit_build(struct bit *tree, int nmemb, double *val);
struct bit *bit_build(int nmemb, double *val);
struct bit *bit_alloc(int nedges);
void bit_update(struct bit *tree, int ix, double diff);
double bit_cumfreq(struct bit *tree, int ix);

static inline double bit_total(struct bit *tree)
{
	return bit_cumfreq(tree, tree->n);
}

double bit_getvalue(struct bit *tree, int ix);
//int bit_getindex_old(struct bit *tree, double prob);
int bit_getindex(struct bit *tree, double freq);
void bit_clear(struct bit *tree);
void bit_free(struct bit *tree);
int bit_append(struct bit *tree, double val);
void bit_print(struct bit *tree);

#endif
