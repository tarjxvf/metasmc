/* Subroutines for binary indexed tree */
#include "bit.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static inline int bit_lowbit_remove(int ix)
{
	return ix & (ix - 1);
}

static inline int bit_lowbit_get(int ix)
{
	return ix & (-ix);
}

void __bit_build_old(struct bit *tree, int nmemb, double *val)
{
	int i, j, maxn;

	maxn = 1;
	while(maxn < nmemb + 1) maxn <<= 1;

	tree->freq = realloc(tree->freq, sizeof(double) * maxn);
	tree->n = nmemb;
	tree->maxnode = maxn;
	tree->freq[0] = 0;
	for(i = 1; i < maxn; i++){
		tree->freq[i] = (i <= nmemb)?val[i - 1]:0;
		for(j = i - 2; j >= i - bit_lowbit_get(i); j--)
			tree->freq[i] += (j < nmemb)?val[j]:0;
	}
}

void __bit_build(struct bit *tree, int nmemb, double *val)
{
	int i, j, maxn, parent;

	maxn = 1;
	while(maxn < nmemb + 1) maxn <<= 1;

	if(maxn > tree->maxnode){
		tree->freq = realloc(tree->freq, sizeof(double) * maxn);
		tree->maxnode = maxn;
	}
	tree->n = nmemb;
	tree->freq[0] = 0;
	/* Forward pass that calculate true cumulative frequency. */
	for(i = 1; i <= nmemb; i++)
		tree->freq[i] = tree->freq[i - 1] + val[i - 1];

//	for(;i < maxn; i++)
//		tree->freq[i] = tree->freq[nmemb];

	/* Backward pass that calculates tree values. */
	for(i = nmemb; i >= 1; i--){
		parent = bit_lowbit_remove(i);
		tree->freq[i] -= tree->freq[parent];
	}
}

struct bit *bit_build(int nmemb, double *val)
{
	struct bit *tree;

	tree = malloc(sizeof(struct bit));
	tree->freq = NULL;
	__bit_build(tree, nmemb, val);

	return tree;
}

/* Allocate an empty binary indexed tree with all frequencies zero. */
struct bit *bit_alloc(int n)
{
	struct bit *tree;
	int maxn;

	maxn = 1;
	while(maxn < n + 1) maxn <<= 1;
	tree = malloc(sizeof(struct bit));
	tree->freq = malloc(sizeof(double) * maxn);
	memset(tree->freq, 0, sizeof(double) * maxn);
	tree->maxnode = maxn;
	tree->n = 0;

	return tree;
}

void bit_update(struct bit *tree, int ix, double diff)
{
	do{
		tree->freq[ix] += diff;
		ix += bit_lowbit_get(ix);
	}while(ix <= tree->n);
}

/* Read cumulative value */
double bit_cumfreq(struct bit *tree, int ix)
{
	double sum;

	sum = tree->freq[0];
	while(ix > 0){
		sum += tree->freq[ix];
		ix = bit_lowbit_remove(ix);
	}

	return sum;
}

/* Read individual value. */
double bit_getvalue(struct bit *tree, int ix)
{
	double val;
	int parent;

	val = tree->freq[ix];
	if(ix > 0){
		parent = bit_lowbit_remove(ix);
		ix--;
		while(ix > parent){
			val -= tree->freq[ix];
			ix = bit_lowbit_remove(ix);
		}
	}

	return val;
}

/* Original algorithm in Fenwick 1994. This is obsoleted because it does not support zero values. */
int bit_getindex(struct bit *tree, double prob)
{
	int index, testix, mask;

	index = 0;
	mask = tree->maxnode >> 1;
	if(prob > tree->freq[0]){
		while(mask > 0){
			testix = index + mask;
			if(testix <= tree->n && prob >= tree->freq[testix]){
				index = testix;
				prob -= tree->freq[index];
			}
			mask >>= 1;
		}
	}

	return index + 1;
}

int bit_getindex_new(struct bit *tree, double freq)
{
	int baseix, testix, half, size;

	if(freq < tree->freq[0])
		return 0;
	baseix = 0;
	freq -= tree->freq[0];
//	freq++;
	size = tree->maxnode;
	half = size >> 1;
	while(half > 0){
		testix = baseix + half;
		if(testix <= tree->n && freq > tree->freq[testix]){
			baseix = testix;
			freq -= tree->freq[baseix];
		}
		half >>= 1;
	}

	return baseix + 1;
}

int bit_append(struct bit *tree, double val)
{
	double oldsum, parent_val;
	int n, parent;

	n = tree->n;
	oldsum = bit_total(tree);
	n++;
	if(n >= tree->maxnode){
//		oldval = bit_total(tree);
		tree->freq = realloc(tree->freq, sizeof(double) * (tree->maxnode << 1));
		memset(tree->freq + tree->maxnode, 0, sizeof(double) * tree->maxnode);
		tree->freq[n] = oldsum + val;
		tree->maxnode <<= 1;

	}else{
		parent_val = bit_cumfreq(tree, bit_lowbit_remove(n));
		tree->freq[n] = oldsum + val - parent_val;
//		bit_update(tree, n, val);
	}

	tree->n = n;

	return n;
}

void bit_clear(struct bit *tree)
{
	memset(tree->freq, 0, sizeof(double) * tree->maxnode);
	tree->n = 0;
}

void bit_free(struct bit *tree)
{
	free(tree->freq);
	free(tree);
}

void bit_print(struct bit *tree)
{
	int i;

	fprintf(stderr, "%s: %d\n", __func__, __LINE__);
	fprintf(stderr, "Tree values=");
	for(i = 1; i <= tree->n; i++)
		fprintf(stderr, "%.6f, ", tree->freq[i]);
	fprintf(stderr, "\n");

	fprintf(stderr, "Cumulative frequencies=");
	for(i = 1; i <= tree->n; i++)
		fprintf(stderr, "%.6f, ", bit_cumfreq(tree, i));
	fprintf(stderr, "\n");

	fprintf(stderr, "Individual value=");
	for(i = 1; i <= tree->n; i++)
		fprintf(stderr, "%.6f, ", bit_getvalue(tree, i));
	fprintf(stderr, "\n");
}

