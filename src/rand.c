#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <unistd.h>
#include "mt19937-64/mt64.h"
//#include "mtwist-1.5/mtwist.h"
//#include "mt19937ar-nrl/mt.h"

FILE *seedfilp;

#define BUF_SIZE 100000
double rnd_buffer[BUF_SIZE];
int buf_cur = BUF_SIZE;

/* Generating [0, 1) uniform random variable. */
/*double dunif01()
{
	int r;
	r = rand();
//	fprintf(seedfilp, "%d\n", r);
	return (double)r / RAND_MAX;
}*/

void populate_buf()
{
	int i;
	for(i = 0; i < BUF_SIZE; i++)
		rnd_buffer[i] = genrand64_real3();
	buf_cur = 0;
}

/* Generate uniform distribution on [0,1). */
double dunif01()
{
	double r;
//	r = mt_drand();
	if(buf_cur >= BUF_SIZE)
		populate_buf();
	r = rnd_buffer[buf_cur++];
//	r = genrand64_real3();

//fprintf(stderr, "unifRV()=%.10f\n", r);
	return r;

//	return genrand_real1();
//	return mt_drand();
//	return genrand64_real1();
}

double dexp(double lambda)
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

unsigned int dunif(unsigned int max)
{
//	unsigned long long r;
//	unsigned long int r;
//	r = genrand_int32();
//	r = mt_lrand();
//	r = genrand64_int64();
//	fprintf(seedfilp, "%d\n", r);
//	return r % max;
	double r = dunif01();
	return r * max;
}

/* This function is from ms. */
double gasdev(double m, double v)
{
	static int iset=0;
	static float gset;
	float fac,r,v1,v2;

	if(iset == 0) {
		do {
			v1=2.0*dunif01()-1.0;
			v2=2.0*dunif01()-1.0;
			r=v1*v1+v2*v2;
		} while (r >= 1.0);
		fac=sqrt(-2.0*log(r)/r);
		gset= v1*fac;
		iset=1;
		return( m + sqrt(v)*v2*fac);
	}else{
		iset=0;
		return( m + sqrt(v)*gset ) ;
	}
}

/* Generating poisson random variable. This function is from ms. */
int poisso(double u)
{
	double  cump, ru, p;
	int i=1;

	if( u > 30. ){
		i =  (int)(0.5 + gasdev(u,u)) ;
		if( i < 0 ) return 0 ;
		else return i ;
	}

	ru = dunif01();
	p = exp(-u);
	if( ru < p) return(0);
	cump = p;

	while( ru > ( cump += (p *= u/i ) ) )
		i++;

	return i;
}

/* Initialize random number generator. */
//void init_rand(unsigned int seed)
void init_rand(unsigned int seed)
{
//	unsigned int seed = 1427717129;
//	seed = 1427211304;
//	seed = time(NULL);
//	init_genrand(seed);
//	mt_seed32new(seed);
	init_genrand64(seed);
//	seedfilp = fopen("rands.txt", "w+");
}

/* Initialize random number generator. */
void seed()
{
//	struct timespec t;
//	clock_gettime(CLOCK_MONOTONIC, &t);
//	mt_seed32new(t.tv_nsec);
//	init_genrand(t.tv_nsec);
//	init_genrand64(t.tv_nsec);
}

void finish_rand()
{
//	fclose(seedfilp);
}

