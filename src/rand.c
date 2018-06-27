#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <unistd.h>

#include "rand.h"
//#include "mtwist-1.5/mtwist.h"
//#include "mt19937ar-nrl/mt.h"

FILE *seedfilp;

dsfmt_t dsfmt;
//uint32_t int_buffer[BUF_SIZE];
double real_buffer[BUF_SIZE];
int buf_cur = 0;

/* Generating [0, 1) uniform random variable. */
/*double dunif01()
{
	int r;
	r = rand();
//	fprintf(seedfilp, "%d\n", r);
	return (double)r / RAND_MAX;
}*/

/*void populate_buf()
{
	int i;
	for(i = 0; i < BUF_SIZE; i++)
		rnd_buffer[i] = genrand64_real3();
	buf_cur = 0;
}*/

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

//	init_genrand64(seed);

	dsfmt_init_gen_rand(&dsfmt, seed);
	dsfmt_fill_array_open_open(&dsfmt, real_buffer, BUF_SIZE);
//	st_fill_array(&dsfmt, int_buffer, BUF_SIZE);

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

