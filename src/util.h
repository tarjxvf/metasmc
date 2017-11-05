#ifndef UTIL_H
#define UTIL_H

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#ifdef __cplusplus
extern "C"{
#endif

#define NEXT_NOBLANK(fp, ch) 	while(isblank((ch) = fgetc(fp)))
#define NEXT_NOBLANK_2(fp, ch)	while(isspace(ch)) ch = fgetc(fp)
#define NEXT_NOSPACE(fp, ch) 	while(isspace((ch) = fgetc(fp)))
#define NEXT_NOSPACE_2(fp, ch) 	while(isspace((ch)) ch = fgetc(fp))

#define MIN(a,b)		((a) <= (b))?(a):(b)
#define MAX(a,b)		((a) >= (b))?(a):(b)

/* Given a list of sorted numbers, search what interval hitted by key. */
void *search_interval(void *pos, void *base, int nmemb, size_t size, int (*compar)(void *, void *));

/* Read an integer from filp. */
int read_integer(FILE *filp, int *val);

/* Read an double from filp. */
int read_double(FILE *filp, double *val);

/* Compare two integers. */
int intcompar(int *a, int *b);

/* Compare two doubles. */
int dblcompar(const double *a, const double *b);

/* base: pointer to the array to be sorted.
 * start: index of the root of the subheap.
 * end: index of final element of the subheap.
 * size: size of each element.
 * data: Additional data for comparison.
 * tmp: working space for swapping elements. */
void heap_restore(void *base, int start, int end, size_t size, int (*compar)(const void *, const void *, const void *), void *data, void *tmp);

/* The parameters are the same as qsort */
void heap_build(void *base, size_t nmemb, size_t size, int (*compar)(const void *, const void *, const void *), void *data);
void heap_sort(void *base, size_t nmemb, size_t size, int (*compar)(const void *, const void *, const void *), void *data);

#ifdef __cplusplus
}
#endif

#endif
