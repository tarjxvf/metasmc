/* This file include mutation generators that modify sequence along genealogy.
 * Supported models include infinite-site model and general reversible model. */
#include <math.h>
#include <string.h>

#include <lapacke.h>

#include "rand.h"
#include "util.h"
#include "mutation.h"

#define LAUNCH_GEEV(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, lcvr, info)			info = LAPACKE_dgeev(LAPACK_COL_MAJOR, (jobvl), (jobvr), (n), (a), (lda), (wr), (wi), (vl), (ldvl), (vr), (lcvr))
#define LAUHCN_GETRF(m, n, a, lda, ipiv, info)							info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, (m), (n), (a), (lda), (ipiv))
#define LAUNCH_GETRI(n, a, lda, ipiv, info)							info = LAPACKE_dgetri(LAPACK_COL_MAJOR, (n), (a), (lda), (ipiv))

char *model_name[] = {"Infinite-site", "JC69", "K80", "F81", "HKY85", "T92", "TN93", "GTR"};
int npar[] = {0, 0, 1, 0, 1, 1, 2, 6};	// Lengths of parameter list for all models

void (*generate_sequence)(struct reference *, struct genealogy *, int, int);

/* Determine whether n2 is descendent of n1. */
int isdesc(struct genealogy *G, struct node *n1, struct node *n2)
{
	while(n1 != n2 && n2 != G->root)
		n2 = n2->in;
	return n1 == n2;
}
/*
int intcompar(const void *a, const void *b)
{
	return *(int *)a - *(int *)b;
}*/

void locate(int from, int to, int segsites, int *pos)
{
	int len, i;

	len = to - from;
	for(i = 0; i < segsites; i++)
		pos[i] = from + dunif(len);
	qsort((void *)pos, (size_t)segsites, sizeof(int), (__compar_fn_t)intcompar);
}

/* Infinite-site model mutation generator. */
void generate_sequence_infinite_fast(struct reference *ref, struct genealogy *G, int from, int to)
{
	struct list_head *el, *fgl;
	struct config *cfg;
	int segsites, pop, curr, i, j, p, r, k, nR, *R;
	double theta, tt, t, *weights;
	struct node **nmut;	/* List of nodes above which mutation event occur */
	struct node *e, **edges;
	char *ances, *deriv;
	int nedges, maxedges;
	struct frag *fgset;
	struct read **rdset, *rds;

#ifdef DEBUG
	fprintf(stderr, "Entering function %s, from=%d, to=%d\n", __func__, from, to);
#endif
	cfg = G->cfg;
	fgset = cfg->prof->fgset;
	rdset = cfg->prof->rdset;

	R = G->R[G->curridx];
	nR = G->nR[G->curridx];

	maxedges = nedges = 0;
	for(i = 0; i < cfg->npop_all; i++)
		maxedges += G->pops[i].nedges;

	// Calculate total rate of mutation by iterating over all edges.
	weights = malloc(sizeof(double) * (maxedges + 1));
	edges = malloc(sizeof(struct edge *) * maxedges);
	tt = 0;
	weights[0] = 0;
//	el = G->e_list.front;
//	while(el){
	for(i = 1; i < G->tr_xover->maxnodes; i++){
		e = G->tr_xover->edges[i];
		if(e){
			double w;

#ifdef DEBUG
			fprintf(stderr, "i=%d, e=%x(%.6f)\n", i, e, e->t);
#endif

//			e = (struct edge *)GET_OBJ(el);
			edges[nedges] = e;
			pop = e->pop;
			w = (e->in->t - e->t) * G->pops[pop].mmut->theta;
			weights[nedges + 1] = weights[nedges] + w;
			nedges++;
//			el = el->next;
		}
	}
	tt = weights[nedges];

	/* Skip characters until from. */
	while(ref->curr < from && !feof(ref->filp)){
		if(nucl_index(fgetc(ref->filp)) >= 0)
			ref->curr++;
	}

	if(feof(ref->filp))
		exit(-1);
#ifdef DEBUG
	fprintf(stderr, "ref->curr=%d, ftell()=%d\n", ref->curr, ftell(ref->filp));
#endif
	ances = malloc(sizeof(char) * (to - from + 1));
	for(i = 0; i < to - from; i++){
		int ch;
		while(nucl_index(ch = fgetc(ref->filp)) < 0);
		ref->curr++;
		ances[i] = ch;
	}
	deriv = malloc(sizeof(char) * (to - from + 1));
	nmut = malloc(sizeof(struct node *) * (to - from));
	memset(nmut, 0, sizeof(struct node *) * (to - from));
	deriv[to - from] = ances[to - from] = '\0';

	/* Generate mutation events. */
	segsites = 0;
	for(i = 0; i < to - from; i++){
		if(dunif01() > exp(-tt)){/* Generate a new mutation. */
			int a, an, k;
			double wsum, u, *hit;

			/* Select point of mutation event. */
//			el = G->e_list.front;
			u = dunif01() * tt;
			// Get mutation point using binary search
			hit = (double *)search_interval(&u, (void *)weights, nedges + 1, sizeof(double), (int (*)(void *, void *))dblcompar);

			nmut[i] = edges[hit - weights];	// First node below mutation event

			/* Generate derived allele. */
			a = dunif(3);
			an = nucl_index(ances[i]);
			for(k = 0; k < 4; k++){
				if(k != an){
					if(a > 0)
						a--;
					else
						break;
				}
			}
			deriv[i] = nucl[k];
			segsites++;

		}else{
			deriv[i] = ances[i];
		}
	}
#ifdef DEBUG
	fprintf(stderr, "%d: ances=%s\n", __LINE__, ances);
	fprintf(stderr, "%d: deriv=%s\n", __LINE__, deriv);
#endif

	for(k = 0; k < G->nR[G->curridx]; k++){
//	fgl = G->r_list->front;
//	while(fgl){
		struct frag *fg;
		struct read *rd;

//		fg = *((struct frag **)GET_OBJ(fgl));
		fg = &fgset[R[k]];
		rds = rdset[R[k]];
		for(r = 0; r < fg->nread; r++){
			int lb, ub, j;

			rd = &rds[r];
			lb = (from < rd->start)?rd->start:from;
			ub = (to > rd->end)?rd->end:to;
			j = lb - from;
#ifdef DEBUG
			fprintf(stderr, "%d: lb=%d, ub=%d, rd->start=%d, rd->end=%d\n", __LINE__, lb, ub, rd->start, rd->end);
#endif
			for(i = 0, p = lb; p < ub && i < rd->end - rd->start; i++, j++, p++){
#ifdef DEBUG
//				fprintf(stderr, "%d: p=%d, j=%d\n", __LINE__, p, j);
#endif
				if(ances[j] != deriv[j]){	/* j is a polymorphic site */
					if(isdesc(G, nmut[j], (struct node *)fg->nd)){
						rd->seq[p - rd->start] = deriv[j];

					}else{
						rd->seq[p - rd->start] = ances[j];
					}

				}else{
					rd->seq[p - rd->start] = ances[j];
				}
			}
		}
//		fgl = fgl->next;
	}

	free(nmut);
	free(deriv);
	free(ances);
	free(weights);
	free(edges);

#ifdef DEBUG
	fprintf(stderr, "%d: ref->curr=%d, ftell()=%d\n", __LINE__, ref->curr, ftell(ref->filp));
	fprintf(stderr, "Finishing %s\n\n", __func__);
#endif
}

void generate_sequence_infinite_slow(struct reference *ref, struct genealogy *G, int from, int to)
{
	struct config *cfg;
	int segsites, *pos, pop, curr, i, j, p, r, nR, *R;
	double theta, tt, t;
	struct node *e;
	char *ances;
	struct frag *fgset;
	struct read **rdset, *rds;

#ifdef DEBUG
	fprintf(stderr, "Entering function %s, from=%d, to=%d\n", __func__, from, to);
#endif
	cfg = G->cfg;
	fgset = cfg->prof->fgset;
	rdset = cfg->prof->rdset;
	R = G->R[G->curridx];
	nR = G->nR[G->curridx];

	theta = G->pops[0].mmut->theta;
	tt = theta * G->total;
	pos = malloc(sizeof(int) * (to - from + 1));
	segsites = 0;
	j = 0;
	for(i = 0; i < to - from; i++){
		int nextpos;
		nextpos = (int)(-log(dunif01()) / tt);
		if(nextpos < to){
			pos[j] = nextpos;
			j++;
		}
	}
	segsites = j;
	pos[j] = to + 1;

#ifdef DEBUG
	fprintf(stderr, "segsites=%d, tt=%.6f, pos=[", segsites, tt);
	for(i = 0; i < segsites; i++)
		fprintf(stderr, "%d, ", pos[i]);
	fprintf(stderr, "]\n", pos[i]);
#endif
	/* Skip characters until from. */
	while(ref->curr < from && !feof(ref->filp)){
		if(nucl_index(fgetc(ref->filp)) >= 0)
			ref->curr++;
	}

	ances = malloc(sizeof(char) * (to - from + 1));
	for(i = 0; i < to - from; i++){
		int ch;
		while(nucl_index(ch = fgetc(ref->filp)) < 0);
		ref->curr++;
		ances[i] = ch;
	}

	p = from;
	j = 0;
	for(i = 0; i < to - from; i++, p++){
		struct list_head *fgl;
		char anc, deriv;

		anc = ances[i];
#ifdef DEBUG
		fprintf(stderr, "%d: position %d, ancestral allele=%c\n", __LINE__, p, anc);
#endif
		if(p == pos[j]){/* Generate mutation. */
			int a, an, k;

			j++;
			/* Determine ancestral allele. */
			a = dunif(3);
			an = nucl_index(anc);
#ifdef DEBUG
			fprintf(stderr, "Generating mutation site, a=%d, an=%d", a, an);
#endif
			for(k = 0; k < 4; k++){
				if(k != an){
					if(a > 0)
						a--;
					else
						break;
				}
			}
			deriv = nucl[k];
#ifdef DEBUG
			fprintf(stderr, ", derived allele=%c\n", deriv);
#endif
			rnd_select_point(G, &e, &pop, &t);	/* Select point of mutation event. */
#ifdef DEBUG
			fprintf(stderr, "Mutation on edge %x in population %d at %.6f\n", e, pop, t);
#endif
//			fgl = G->r_list->front;
//			while(fgl){
			for(i = 0; i < nR; i++){
				struct frag *fg;
				struct read *rd;

//				fg = *((struct frag **)GET_OBJ(fgl));
				fg = &fgset[R[i]];
				rds = rdset[R[i]];
				for(r = 0; r < fg->nread; r++){
					rd = &rds[r];
					if(rd->start <= p && p < rd->end){
						/* Generate derived allele at random. */
						if(isdesc(G, e, (struct node *)fg->nd)){
#ifdef DEBUG
							fprintf(stderr, "%x is descendent of %x\n", fg->nd, e);
#endif
							rd->seq[p - rd->start] = deriv;

						}else{
#ifdef DEBUG
							fprintf(stderr, "%x is not descendent of %x\n", fg->nd, e);
#endif
							rd->seq[p - rd->start] = anc;
						}
					}
				}
//				fgl = fgl->next;
			}

		}else{	/* Copy reference allele. */
//			fgl = G->r_list->front;
//			while(fgl){
			for(i = 0; i < nR; i++){
				struct frag *fg;
				struct read *rd;

//				fg = *((struct frag **)GET_OBJ(fgl));
				fg = &fgset[R[i]];
				rds = rdset[R[i]];
#ifdef DEBUG
				fprintf(stderr, "Processing fragment %d(node %x), nreads=%d\n", fg->id, fg->nd, fg->nread);
#endif
				for(r = 0; r < fg->nread; r++){
					rd = &rds[r];
					if(rd->start <= p && p < rd->end){
						rd->seq = rd->seq;
						rd->seq[p - rd->start] = anc;
					}
				}
				fgl = fgl->next;
			}
		}
	}

	free(ances);
	free(pos);
#ifdef DEBUG
	fprintf(stderr, "Finishing %s\n\n", __func__);
#endif
}

/* Eigenvalue decomposition.
 * A: (input) nxn matrix.
 * Q: (output) nxn matrix of which each column is an eigenvector.
 * D: (output) A n-dimensional array of eigenvalues. 
 * Qinv: (output) inverse of Q */
void eigen_decompose(int n, double *A, double *Q, double *D, double *Qinv)
{
	double *wi, *vl;
	int *ipiv, info, i;
	ipiv = (int *)malloc(sizeof(int) * n);

	wi = (double *)malloc(sizeof(double) * n);

        LAUNCH_GEEV('N', 'V', n, A, n, D, wi, NULL, n, Q, n, info);

	memcpy(Qinv, Q, sizeof(double) * n * n);
	LAUHCN_GETRF(n, n, Qinv, n, ipiv, info);
	LAUNCH_GETRI(n, Qinv, n, ipiv, info);
	free(wi);
	free(ipiv);
}

/* Set Q matrix of mutation models and compute its eigendecomposition. */
void init_mutation_model(struct mutation *mmut)
{
	double Q[SQNUCS], P[SQNUCS], Pinv[SQNUCS];// Mutation matrix Q
	int i, j, k;

#ifdef DEBUG
	fprintf(stderr, "Setting model ");
#endif
	/* Setup Q matrix. */
	memset(Q, 0, sizeof(double) * SQNUCS);
	if(mmut->model == MODEL_JC69){
		for(i = 0; i < NUM_NUCS; i++)
			for(j = i + 1; j < NUM_NUCS; j++)
				Q[j * NUM_NUCS + i] = Q[i * NUM_NUCS + j] = 0.25;

	}else if(mmut->model == MODEL_K80){
		Q[0 * NUM_NUCS + 2] = Q[2 * NUM_NUCS + 0] = Q[1 * NUM_NUCS + 3] = Q[3 * NUM_NUCS + 1] = mmut->mpar[0];
		Q[0 * NUM_NUCS + 1] = Q[0 * NUM_NUCS + 3] = Q[1 * NUM_NUCS + 0] = Q[1 * NUM_NUCS + 2] = Q[2 * NUM_NUCS + 1] = Q[2 * NUM_NUCS + 3] = Q[3 * NUM_NUCS + 0] = Q[3 * NUM_NUCS + 2] = 1;

	}else if(mmut->model == MODEL_F81){
		for(i = 0; i < NUM_NUCS; i++)
			for(j = 0; j < NUM_NUCS; j++)
				if(i != j)
					Q[i * NUM_NUCS + j] = mmut->pi[j];

	}else if(mmut->model == MODEL_HKY85){
		Q[0 * NUM_NUCS + 2] = Q[2 * NUM_NUCS + 0] = Q[1 * NUM_NUCS + 3] = Q[3 * NUM_NUCS + 1] = mmut->mpar[0];
		Q[0 * NUM_NUCS + 1] = Q[0 * NUM_NUCS + 3] = Q[1 * NUM_NUCS + 0] = Q[1 * NUM_NUCS + 2] = Q[2 * NUM_NUCS + 1] = Q[2 * NUM_NUCS + 3] = Q[3 * NUM_NUCS + 0] = Q[3 * NUM_NUCS + 2] = 1;

		for(i = 0; i < NUM_NUCS; i++)
			for(j = 0; j < NUM_NUCS; j++)
				Q[i * NUM_NUCS + j] *= mmut->pi[j];

	}else if(mmut->model == MODEL_T92){
		for(i = 0; i < NUM_NUCS; i++){
			for(j = 0; j < NUM_NUCS; j++){
				if(i != j){
				}
			}
		}

	}else if(mmut->model == MODEL_TN93){
		Q[0 * NUM_NUCS + 2] = Q[2 * NUM_NUCS + 0] = mmut->mpar[0];
		Q[1 * NUM_NUCS + 3] = Q[3 * NUM_NUCS + 1] = mmut->mpar[1];
		Q[0 * NUM_NUCS + 1] = Q[0 * NUM_NUCS + 3] = Q[1 * NUM_NUCS + 0] = Q[1 * NUM_NUCS + 2] = Q[2 * NUM_NUCS + 1] = Q[2 * NUM_NUCS + 3] = Q[3 * NUM_NUCS + 0] = Q[3 * NUM_NUCS + 2] = 1;

		for(i = 0; i < NUM_NUCS; i++)
			for(j = 0; j < NUM_NUCS; j++)
				Q[i * NUM_NUCS + j] *= mmut->pi[j];

	}else if(mmut->model == MODEL_GTR){
		int k;

#ifdef DEBUG
		fprintf(stderr, "GTR\n");
#endif

		k = 0;
		for(i = 0; i < NUM_NUCS; i++)
			for(j = i + 1; j < NUM_NUCS; j++){
				Q[j * NUM_NUCS + i] = Q[i * NUM_NUCS + j] = mmut->mpar[k];
				k++;
			}

		for(i = 0; i < NUM_NUCS; i++)
			for(j = 0; j < NUM_NUCS; j++)
				Q[i * NUM_NUCS + j] *= mmut->pi[j];
	}

	for(i = 0; i < NUM_NUCS; i++)
		for(j = 0; j < NUM_NUCS; j++)
			if(i != j)
				Q[i * NUM_NUCS + i] -= Q[i * NUM_NUCS + j];
#ifdef DEBUG
	fprintf(stderr, "Q matrix:[\n");
	for(i = 0; i < NUM_NUCS; i++){
		for(j = 0; j < NUM_NUCS; j++)
			fprintf(stderr, "%.6f, ", Q[i * NUM_NUCS + j]);
		fprintf(stderr, "\n");
	}
	fprintf(stderr, "\n");
#endif
	eigen_decompose(NUM_NUCS, Q, P, mmut->D, Pinv);
	for(i = 0; i < NUM_NUCS; i++)
		for(j = 0; j < NUM_NUCS; j++)
			for(k = 0; k < NUM_NUCS; k++)
				mmut->Cijk[k][j][i] = P[j * NUM_NUCS + i] * Pinv[k * NUM_NUCS + j];	// Note that P matrix is column-major
#ifdef DEBUG
	fprintf(stderr, "Cik:[\n");
	for(i = 0; i < NUM_NUCS; i++){
		for(k = 0; k < NUM_NUCS; k++){
			fprintf(stderr, "[", j);
			for(j = 0; j < NUM_NUCS; j++){
				fprintf(stderr, "%.6f, ", mmut->Cijk[k][j][i]);
			}
			fprintf(stderr, "]\n");
		}
	}
	fprintf(stderr, "]\n");

	fprintf(stderr, "Eigenvalues of Q matrix:[");
	for(i = 0; i < NUM_NUCS; i++)
		fprintf(stderr, "%.6f, ", mmut->D[i]);
	fprintf(stderr, "]\n");

	fprintf(stderr, "P matrix:[\n");
	for(i = 0; i < NUM_NUCS; i++){
		for(j = 0; j < NUM_NUCS; j++)
			fprintf(stderr, "%.6f, ", P[i * NUM_NUCS + j]);
		fprintf(stderr, "\n");
	}
	fprintf(stderr, "]\n");
	fprintf(stderr, "Inverse of P matrix:[\n");
	for(i = 0; i < NUM_NUCS; i++){
		for(j = 0; j < NUM_NUCS; j++)
			fprintf(stderr, "%.6f, ", Pinv[i * NUM_NUCS + j]);
		fprintf(stderr, "\n");
	}
	fprintf(stderr, "]\n");
#endif
}

void dump_mutation_model(struct mutation *mmut)
{
	int i;

	fprintf(stderr, " Model: %s, mutation rate=%.6f\n  mpar: [", model_name[mmut->model], mmut->theta);
	for(i = 0; i < npar[mmut->model]; i++){
		fprintf(stderr, "%.6f, ", mmut->mpar[i]);
	}
	fprintf(stderr, "]\n");
	fprintf(stderr, "  pi:[%.6f, %.6f, %.6f, %.6f]\n", mmut->pi[0], mmut->pi[1], mmut->pi[2], mmut->pi[3]);
}

void setmatrix(struct mutation *mmut, double *Pmut, double t)
{
	double beta, ebeta, expD[NUM_NUCS];
	int i, j, k;

#ifdef DEBUG
	fprintf(stderr, "Setting mutation probabilities, t=%.6f, theta=%.6f\n", t, mmut->theta);
#endif

	if(mmut->theta == 0 || t == 0){
		memset(Pmut, 0, sizeof(double) * SQNUCS);
		for(i = 0; i < NUM_NUCS; i++)
			Pmut[i * NUM_NUCS + i] = 1;
		return;
	}

	/* Compute close-form solutions if available. */
	switch(mmut->model){
		case 0:	// Infinite-site model
			return;

		case MODEL_JC69:
#ifdef DEBUG
			fprintf(stderr, "Close-form solution of JC69 model:[\n");
			for(i = 0; i < NUM_NUCS; i++){
				for(j = 0; j < NUM_NUCS; j++)
					if(i == j)
						fprintf(stderr, "%.6f, ", 0.25 + 0.75 * exp(-mmut->theta * t));
					else
						fprintf(stderr, "%.6f, ", 0.25 - 0.25 * exp(-mmut->theta * t));
				fprintf(stderr, "\n");
			}
			fprintf(stderr, "]\n");
#endif
			break;

		case MODEL_K80:
			break;

		case MODEL_F81:
			break;

		case MODEL_HKY85:
			break;

		case MODEL_T92:
			break;

		case MODEL_TN93:
			break;

		case MODEL_GTR:
			break;
	}

	for(j = 0; j < NUM_NUCS; j++)
		expD[j] = exp(mmut->D[j] * mmut->theta * t);

	memset(Pmut, 0, sizeof(double) * SQNUCS);
	for(i = 0; i < NUM_NUCS; i++)
		for(k = 0; k < NUM_NUCS; k++)
			for(j = 0; j < NUM_NUCS; j++)
				Pmut[i * NUM_NUCS + k] += mmut->Cijk[i][j][k] * expD[j];
#ifdef DEBUG
	fprintf(stderr, "Mutation probabilities:[\n");
	for(i = 0; i < NUM_NUCS; i++){
		for(j = 0; j < NUM_NUCS; j++)
			fprintf(stderr, "%.6f, ", Pmut[i * NUM_NUCS + j]);
		fprintf(stderr, "\n");
	}
	fprintf(stderr, "]\n");
#endif
}

void mutate(struct mutation *mmut, double *Pmut, char *seq, int len)
{
	double u;
	int i, j;
#ifdef DEBUG
	fprintf(stderr, "Entering function %s\n", __func__);
	fprintf(stderr, "Mutation probabilities:[\n");
	for(i = 0; i < NUM_NUCS; i++){
		for(j = 0; j < NUM_NUCS; j++)
			fprintf(stderr, "%.6f, ", Pmut[i * NUM_NUCS + j]);
		fprintf(stderr, "\n");
	}
	fprintf(stderr, "]\n");
#endif
	/* Generate mutated sequence. */
	for(i = 0; i < len; i++){
		int anc;

		anc = nucl_index(seq[i]);
		u = dunif01();
		for(j = 0; j < NUM_NUCS; j++){
			if(u < Pmut[anc * NUM_NUCS + j])
				break;
			u -= Pmut[anc * NUM_NUCS + j];
		}

		seq[i] = nucl[j];
	}
#ifdef DEBUG
	fprintf(stderr, "Finishing function %s\n\n", __func__);
#endif
}

void genseq_model(struct genealogy *G, struct mutation *mmut, struct node *n, char *seq, int from, int to);

void __genseq(struct genealogy *G, struct mutation *mmut, struct node *n, char *seq_old, struct node *n2, int from, int to)
{
	double Pmut[SQNUCS];
	char *seq;
	double dt;

	while(ismigrnode(n2))
		n2 = ((struct migr_node *)n2)->out;

	seq = malloc(sizeof(char) * (to - from));
	memcpy(seq, seq_old, sizeof(char) * (to - from));
	dt = n->t - n2->t;
	setmatrix(mmut, Pmut, dt);
	mutate(mmut, Pmut, seq, to - from);

	genseq_model(G, mmut, n2, seq, from, to);
	free(seq);
}

void genseq_model(struct genealogy *G, struct mutation *mmut, struct node *n, char *seq_old, int from, int to)
{
	struct node *n2;

#ifdef DEBUG
	fprintf(stderr, "Entering %s, n=%x, n->type=%d, from=%d, to=%d\n", __func__, n, n->type, from, to);
#endif
	if(iscoalnode(n)){
		__genseq(G, mmut, n, seq_old, ((struct coal_node *)n)->out[0], from, to);
		__genseq(G, mmut, n, seq_old, ((struct coal_node *)n)->out[1], from, to);

	}else if(issamnode(n)){
		struct sam_node *nd;
		struct frag *fg;
		struct read *rd, *rds, **rdset;
		int i, j, k, f;

		nd = (struct sam_node *)n;
		f = nd->fgid;
		fg = &G->cfg->prof->fgset[f];
		rds = G->cfg->prof->rdset[f];
#ifdef DEBUG
		fprintf(stderr, "frag id=%d\n", nd->fg->id);
#endif
		for(i = 0; i < fg->nread; i++){
			int lb, ub;
//			rd = &fg->rd[i];
			rd = &rds[i];
			lb = (from < rd->start)?rd->start:from;
			ub = (to > rd->end)?rd->end:to;
			for(j = lb; j < ub; j++)
				rd->seq[j - rd->start] = seq_old[j - from];
		}
	}
}

void generate_sequence_model_slow(struct reference *ref, struct genealogy *G, int from, int to)
{
	struct config *cfg;
	double Pmut[SQNUCS];
	char *seq;
	int i;

	cfg = G->cfg;
	/* Skip characters until from. */
	while(ref->curr < from && !feof(ref->filp)){
		if(nucl_index(fgetc(ref->filp)) >= 0)
			ref->curr++;
	}

	if(feof(ref->filp))
		exit(-1);
#ifdef DEBUG
	fprintf(stderr, "ref->curr=%d, ftell()=%d, from=%d, to=%d\n", ref->curr, ftell(ref->filp), from, to);
#endif
	seq = malloc(sizeof(char) * (to - from));
	for(i = 0; i < to - from; i++){
		int ch;
		while(nucl_index(ch = fgetc(ref->filp)) < 0);
		ref->curr++;
		seq[i] = ch;
	}

	if(cfg->tdiv > 0){/* Generate mutated root sequence. */
		double dt;

		dt = cfg->tdiv;
		setmatrix(G->pops[0].mmut, Pmut, dt);
		mutate(G->pops[0].mmut, Pmut, seq, to - from);
	}
	genseq_model(G, G->pops[0].mmut, G->root, seq, from, to);
	free(seq);
#ifdef DEBUG
	fprintf(stderr, "Finishing function %s\n", __func__);
#endif
}

// Check whether the sequence of a coalescent node has already been generated in [off_lb, off_ub)
int mapped(struct coal_node *n, int off_lb, int off_ub)
{
	int m1, m2, ncell, i, j;
	size_t mask;

	if(n->mapped == NULL)
		return 0;

	m1 = off_lb / CELLSIZE;
	m2 = off_ub / CELLSIZE;

	mask = 0x80000000;
	for(i = NCELL_PER_MAP; i > m1 % NCELL_PER_MAP; i--){
		mask >>= 1;
		mask |= 0x80000000;
	}
	if(!(n->mapped[m1 / NCELL_PER_MAP] & mask))
		return 0;

	mask = 0x1;
	for(i = 0; i < m2 % NCELL_PER_MAP; i++){
		mask <<= 1;
		mask |= 0x1;
	}
	if(!(n->mapped[m2 / NCELL_PER_MAP] & mask))
		return 0;

	mask = 0xFFFFFFFF;
	for(i = m1 + 1; i < m2; i++)
		if(!(n->mapped[i / NCELL_PER_MAP] & mask))
			return 0;

	return 1;
}

// Evolve coalescent node
void evolve(struct mutation *mmut, struct coal_node *n2, struct coal_node *n1, int off_lb, int off_ub)
{
	double dt, Pmut[SQNUCS];
	int m1, m2, ncell, i, j, k;
	size_t mask;

	m1 = off_lb / CELLSIZE;
	m2 = off_ub / CELLSIZE;
#ifdef DEBUG
	fprintf(stderr, "Entering function %s: n2=%x, n1=%x, off_lb=%d, m1=%d, off_ub=%d, m2=%d\n", __func__, n2, n1, off_lb, m1, off_ub, m2);
#endif
	dt = n2->t - n1->t;
	setmatrix(mmut, Pmut, dt);
	mask = 1;
	for(j = 0; j < m1 % NCELL_PER_MAP; j++)
		mask <<= 1;
	j = m1 / NCELL_PER_MAP;
	for(i = m1; i <= m2; i++){
#ifdef DEBUG
		fprintf(stderr, "%d: i=%d, j=%d, mapped[j]=%d, mask=%u\n", __LINE__, i, j, n1->mapped[j], mask);
#endif
		if(mask == 0){
			mask = 1;
			j++;
		}
		if(!(n1->mapped[j] & mask)){
			char *seq1, *seq2;

			seq1 = n1->seq + CELLSIZE * i;
			seq2 = n2->seq + CELLSIZE * i;
			memcpy(seq1, seq2, sizeof(char) * CELLSIZE);
			mutate(mmut, Pmut, seq1, CELLSIZE);
			n1->mapped[j] |= mask;
		}
		mask <<= 1;
	}
}

void clean_mapped(struct node *n)
{
	if(iscoalnode(n)){
		struct coal_node *nc;

		nc = AS_COAL_NODE(n);
		clean_mapped(nc->out[0]);
		clean_mapped(nc->out[1]);
		free(nc->seq);
		free(nc->mapped);
		nc->seq = NULL;
		nc->mapped = NULL;

	}else if(issamnode(n)){
		return;

	}else if(ismigrnode(n)){
		struct node *np;

		np = n;
		while(ismigrnode(np)) np = AS_MIGR_NODE(np)->out;
		clean_mapped(np);
	}
}

void generate_sequence_model_fast(struct reference *ref, struct genealogy *G, int from, int to)
{
	struct config *cfg;
	int pop, ncell, nmap, lb, ub, i, j, p, r, nR, *R;
	double dt, Pmut[SQNUCS];
	struct read *rd, *rds;
	char *seq;
	struct frag *fgset;
	struct read **rdset;

	cfg = G->cfg;
	fgset = cfg->prof->fgset;
	rdset = cfg->prof->rdset;

	R = G->R[G->curridx];
	nR = G->nR[G->curridx];

	/* Skip characters until from. */
	while(ref->curr < from && !feof(ref->filp)){
		if(nucl_index(fgetc(ref->filp)) >= 0)
			ref->curr++;
	}

	ncell = (to - from + CELLSIZE - 1) / CELLSIZE;
	nmap = (ncell + NCELL_PER_MAP - 1) / NCELL_PER_MAP;
#ifdef DEBUG
	fprintf(stderr, "Entering function %s: from=%d, to=%d, ncell=%d\n", __func__, from, to, ncell);
#endif
	seq = malloc(sizeof(char) * ncell * CELLSIZE);
	memset(seq, 0, sizeof(char) * ncell * CELLSIZE);
	for(i = 0; i < to - from; i++){
		int ch;
		while(nucl_index(ch = fgetc(ref->filp)) < 0);
		ref->curr++;
		seq[i] = ch;
	}

//	if(cfg->tdiv > 0){/* Generate mutated root sequence. */
//		dt = cfg->tdiv;
//		setmatrix(cfg->mmut[0], Pmut, dt);
//		mutate(cfg->mmut[0], Pmut, seq, to - from);
//	}

	if(iscoalnode(G->root)){
		struct list_head *fgl;
		struct list stack;
#ifdef DEBUG
		fprintf(stderr, "%d: root=%x\n", __LINE__, G->root);
#endif
		((struct coal_node *)G->root)->seq = seq;
		((struct coal_node *)G->root)->mapped = malloc(sizeof(map_t) * nmap);
//		((struct coal_node *)G->root)->mapped = malloc(sizeof(map_t) * ncell);
		for(i = 0; i < nmap; i++)
			((struct coal_node *)G->root)->mapped[i] = 0xFFFFFFFF;

		list_init(&stack);
//		fgl = G->r_list->front;
		// Trace from leaf node upward until a root or coalescent node which is already mapped
//		while(fgl){
		for(i = 0; i < nR; i++){
			struct frag *fg;
			struct read *rd;

//			fg = *((struct frag **)GET_OBJ(fgl));
			fg = &fgset[R[i]];
			rds = rdset[R[i]];
			for(r = 0; r < fg->nread; r++){	// Generate sequence for each read covering current region
				int lb, ub, j;

				rd = &rds[r];
				lb = (from < rd->start)?rd->start:from;	// Start position is either segment left boundary or read start position
				ub = (to > rd->end)?rd->end:to;	// End position is either segment right boundary or read end position

				if(lb < ub){	// This node has to be simulated because its sequences covers current range
					struct coal_node *n2;
					struct node *n1;
					int off_ub, off_lb;

					off_lb = lb - from;
					off_ub = ub - from;
					n1 = (struct node *)fg->nd;
					n2 = (struct coal_node *)fg->nd->in;

					// Find nearest coalescent node (the node right above leaf may be migration node)
					while((struct node *)n2 != G->root && !iscoalnode((struct node *)n2)) n2 = (struct coal_node *)n2->in;

					while(!mapped(n2, off_lb, off_ub)){	// If the sequence of current coalescent has been generated in the region of interest, then stop tracing. Otherwise, push current node to stack and trace upward.
						struct list_head *top;

						// Push current coalescent node (n2) into the stack
						if(n2->seq == NULL){
							n2->seq = (char *)malloc(sizeof(char) * ncell * CELLSIZE);
							memset(n2->seq, 0, sizeof(char) * ncell * CELLSIZE);
						}

						if(n2->mapped == NULL){
//							n2->mapped = malloc(sizeof(map_t) * ncell);
							n2->mapped = malloc(sizeof(map_t) * nmap);
							memset(n2->mapped, 0, sizeof(map_t) * nmap);
						}

						top = malloc(sizeof(struct list_head) + sizeof(struct coal_node *));
						*((struct coal_node **)GET_OBJ(top)) = n2;
						list_add(&stack, GET_OBJ(top));
						n1 = (struct node *)n2;
						n2 = (struct coal_node *)n1->in;

						// Find next coalescent node
						while((struct node *)n2 != G->root && !iscoalnode((struct node *)n2)) n2 = (struct coal_node *)n2->in;
					}

					while(stack.front){
						struct list_head *top;

						// Get top element in the stack
						top = stack.front;
						__list_remove(&stack, top);
						n1 = *((struct node **)GET_OBJ(top));
						evolve(G->pops[n1->pop].mmut, n2, (struct coal_node *)n1, off_lb, off_ub);
						n2 = (struct coal_node *)n1;
						free(top);
					}

					// Evolve leaf node
					dt = n2->t - fg->nd->t;
					setmatrix(G->pops[fg->nd->pop].mmut, Pmut, dt);
#ifdef DEBUG
					fprintf(stderr, "%d: n2=%x, n=%x, rd->start=%d, lb=%d, rd->end=%d, ub=%d\n", __LINE__, n2, n, rd->start, lb, rd->end, ub);
#endif
					for(j = lb; j < ub; j++)
						rd->seq[j - rd->start] = n2->seq[j - from];
//					fprintf(stderr, "%d: n=0x%x\n", __LINE__, n);
					mutate(G->pops[fg->nd->pop].mmut, Pmut, rd->seq + (lb - rd->start), ub - lb);
				}
			}
//			fgl = fgl->next;
		}

	}else if(issamnode(G->root)){	// If root is a sample node, do the following operation (theoretically, you shouldn't go here)
		struct sam_node *n;
		struct frag *fg;

		n = (struct sam_node *)G->root;
		fg = &fgset[n->fgid];
		rds = rdset[n->fgid];
		for(r = 0; r < fg->nread; r++){
			rd = &rds[r];
			rd = &rds[r];
			lb = (from < rd->start)?rd->start:from;
			ub = (to > rd->end)?rd->end:to;
			if(lb < ub){
				for(j = lb; j < ub; j++)
					rd->seq[j - rd->start] = seq[j - from];
			}
		}
		free(seq);
	}

	/* Free mapped and seq at each coalescent nodes. */
	clean_mapped(G->root);
}

