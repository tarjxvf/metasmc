#include "util.h"

#include <string.h>

/* Search interval by the position. */
void *search_interval(void *pos, void *base, int nmemb, size_t size, int (*compar)(void *, void *))
{
	int start, end, pivot, ret;

//fprintf(stderr, "Entering %s: base=%x, nmemb=%d\n", __func__, base, nmemb);
	start = 0;
	end = nmemb;
	pivot = (start + end + 1) / 2;
	while(start != end && pivot < nmemb){
		if((ret = compar(pos, base + pivot * size)) < 0){
			end = pivot - 1;

		}else if(ret > 0){
			start = pivot;

		}else{
			return base + pivot * size;
		}
		pivot = (start + end + 1) / 2;
	}

//fprintf(stderr, "%d: ptr=%x, start=%d, pivot=%d, end=%d\n", __LINE__, base + start * size, start, pivot, end);
	return base + start * size;
//	return base + pivot * size;
}

int read_integer(FILE *filp, int *val)
{
	int ch;

	NEXT_NOBLANK(filp, ch);
	if(!isdigit(ch)){
		fprintf(stderr, "Syntax error: invalid integer.\n");
		return -1;
	}
	*val = (ch - '0');
	while(isdigit(ch = fgetc(filp))) *val = *val * 10 + (ch - '0');
	return ch;
}

int read_double(FILE *filp, double *val)
{
	int ch;

	NEXT_NOBLANK(filp, ch);
	if(!isdigit(ch)){
		fprintf(stderr, "Syntax error: invalid real value.\n");
		return -1;
	}

	*val = ch - '0';
	while(isdigit(ch = fgetc(filp))) *val = *val * 10 + (ch - '0');
	if(ch == '.'){
		// There is fractional part, read it
		double frac = 1. / 10;
		while(isdigit(ch = fgetc(filp))){
			*val += frac * (ch - '0');
			frac /= 10;
		}
	}

	if(ch == 'e' || ch == 'E'){
		int sign, expo, i;

		sign = fgetc(filp);
		ch = read_integer(filp, &expo);

		if(ch < 0)
			return -1;

		if(sign == '-'){
			for(i = 0; i < expo; i++)
				*val /= 10;

		}else if(sign == '+'){
			for(i = 0; i < expo; i++)
				*val *= 10;

		}else{
			fprintf(stderr, "Expect for '+' or '-'\n");
			return -1;
		}
	}

	return ch;
}

/* Compare two integers. */
int intcompar(int *a, int *b)
{
	return *a - *b;
}

int dblcompar(const double *a, const double *b)
{
	if(*a < *b) return -1;
	else if(*a == *b) return 0;
	else return 1;
}

void heap_restore(void *base, int start, int end, size_t size, int (*compar)(const void *, const void *, const void *), void *data, void *tmp)
{
	int j, left, right, imax;

	j = start;
	while(j * 2 + 1 <= end){	// adjust heap
		left = j * 2 + 1;
		right = left + 1;
		imax = j;
		if(compar(base + size * left, base + size * j, data) > 0)
			imax = left;

		if(right <= end && compar(base + size * right, base + size * imax, data) > 0)
			imax = right;

		if(imax == j){
			break;

		}else{
			memcpy(tmp, base + size * j, size);
			memcpy(base + size * j, base + size * imax, size);
			memcpy(base + size * imax, tmp, size);
			j = imax;
		}
	}
}

void heap_build(void *base, size_t nmemb, size_t size, int (*compar)(const void *, const void *, const void *), void *data)
{
	int i, j, par;
	void *tmp;

	tmp = malloc(size);
	for(i = (nmemb - 2) / 2; i >= 0; i--)
		heap_restore(base, i, nmemb - 1, size, compar, data, tmp);
	free(tmp);
}

void heap_sort(void *base, size_t nmemb, size_t size, int (*compar)(const void *, const void *, const void *), void *data)
{
	int i, j, left, right, imax;
	void *tmp;

	tmp = malloc(size);
	heap_build(base, nmemb, size, compar, data);

	for(i = 1; i < nmemb; i++){
		memcpy(tmp, base + size * (nmemb - i), size);
		memcpy(base + size * (nmemb - i), base, size);
		memcpy(base, tmp, size);
		heap_restore(base, 0, nmemb - i - 1, size, compar, data, tmp);
	}
	free(tmp);
}

