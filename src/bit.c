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

struct bit *bit_build(int nmemb, size_t size, double *val)
{
	struct bit *tree;
	int i, j, maxn;

	maxn = 1;
	while(maxn < nmemb + 1) maxn <<= 1;

	tree = malloc(sizeof(struct bit));
	tree->freq = malloc(sizeof(double) * maxn);
	tree->n = nmemb;
	tree->maxnode = maxn;
	tree->freq[0] = 0;
	for(i = 1; i < maxn; i++){
		tree->freq[i] = (i <= nmemb)?val[i - 1]:0;
		for(j = i - 2; j >= i - bit_lowbit_get(i); j--)
			tree->freq[i] += (j < nmemb)?val[j]:0;
	}

	return tree;
}

void bit_update(struct bit *tree, int ix, double diff)
{
	do{
		tree->freq[ix] += diff;
		ix += bit_lowbit_get(ix);
	}while(ix < tree->maxnode);
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
/*int bit_getindex_old(struct bit *tree, double prob)
{
	int index, testix, mask;

	index = 0;
	mask = tree->maxnode >> 1;
//	while(mask > tree->n) mask >>= 1;
	if(prob > tree->freq[0]){
		while(mask > 0){
			testix = index + mask;
			if(prob >= tree->freq[testix]){
				index = testix;
				prob -= tree->freq[index];
			}
			mask >>= 1;
		}
	}

	return index;
}*/

int bit_getindex(struct bit *tree, double freq)
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
		if(freq > tree->freq[testix]){
			baseix = testix;
			freq -= tree->freq[baseix];
		}
		half >>= 1;
	}

	return baseix + 1;
}

void bit_append(struct bit *tree, double val)
{
	double oldval;
	int n, parent;

	n = tree->n;
	n++;
	if(n >= tree->maxnode){
		double oldval;

		oldval = bit_total(tree);
		tree->freq = realloc(tree->freq, sizeof(double) * (tree->maxnode << 1));
		memset(tree->freq + tree->maxnode, 0, sizeof(double) * tree->maxnode);
		tree->freq[n] = oldval + val;
		tree->maxnode <<= 1;

	}else{
		bit_update(tree, n, val);
	}

	tree->n = n;
}

void bit_free(struct bit *tree)
{
	free(tree->freq);
	free(tree);
}

