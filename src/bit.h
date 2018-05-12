#ifndef BIT_H
#define BIT_H

#include <stdio.h>

struct bit {
	int maxnode;
	int n;
	double *freq;
};

struct bit *bit_build(int nmemb, size_t size, double *val);
void bit_update(struct bit *tree, int ix, double diff);
double bit_cumfreq(struct bit *tree, int ix);

static inline double bit_total(struct bit *tree)
{
	return bit_cumfreq(tree, tree->n);
}

double bit_getvalue(struct bit *tree, int ix);
//int bit_getindex_old(struct bit *tree, double prob);
int bit_getindex(struct bit *tree, double freq);
void bit_free(struct bit *tree);
void bit_append(struct bit *tree, double val);

#endif
