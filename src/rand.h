#ifndef RAND_H
#define RAND_H

#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <inttypes.h>
#include "dSFMT.h"

//#include "mt19937-64/mt64.h"

#define BUF_SIZE 8192

extern dsfmt_t dsfmt;
//extern uint32_t int_buffer[BUF_SIZE];
extern double real_buffer[BUF_SIZE];
extern int buf_cur;

/* Generate uniform distribution on [0,1). */
static inline double dunif01()
{
	double r;
//	r = mt_drand();
	if(buf_cur >= BUF_SIZE){
		dsfmt_fill_array_open_open(&dsfmt, real_buffer, BUF_SIZE);
		buf_cur = 0;
	}
	r = real_buffer[buf_cur++];
//	r = genrand64_real3();

//fprintf(stderr, "unifRV()=%.10f\n", r);
	return r;

//	return genrand_real1();
//	return mt_drand();
//	return genrand64_real1();
}

static inline double dexp(double lambda)
{
	return -log(1 - dunif01()) / lambda;
}

/* Generating discrete uniform random variable ranging from 0 to max - 1. */
/*int dunif(int max)
{
	int r;
	r = rand();
	fprintf(seedfilp, "%d\n", r);
	return r % max;
}*/

static inline unsigned int dunif(unsigned int max)
{
	uint32_t r;

//	unsigned long long r;
//	unsigned long int r;
//	r = genrand_int32();
//	r = mt_lrand();
//	r = genrand64_int64();
//	fprintf(seedfilp, "%d\n", r);
//	return r % max;
//	double r = dunif01();
//	if(buf_cur >= BUF_SIZE)
//		sfmt_fill_array32(&sfmt, int_buffer, BUF_SIZE);
//	r = int_buffer[buf_cur++];

//	return r % max;
//	return dunif01() * max;
	r = dsfmt_genrand_uint32(&dsfmt);
	return r % max;
}

void seed();
void finish_rand();
void init_rand(unsigned int seed);

#endif
