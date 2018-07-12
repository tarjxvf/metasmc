#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include "rand.h"
#include "util.h"
#include "global.h"
#include "smc.h"
#include "mutation.h"
#include "evindex.h"

size_t evsize[] = {sizeof(struct coal_event), sizeof(struct migr_event), sizeof(struct grow_event), sizeof(struct size_event), sizeof(struct rmig_event), sizeof(struct gmig_event), sizeof(struct gsiz_event), sizeof(struct ggro_event), sizeof(struct join_event), sizeof(struct splt_event), sizeof(struct event), sizeof(struct event), sizeof(struct samp_event)};

//struct event *alloc_event(struct config *cf, int type, int pop, double t)
struct event *alloc_event(struct config *cfg, int type, double t)
{
	struct list_head *l;
	struct event *ev;
	int npop_all;

	npop_all = cfg->npop + cfg->nsplt;
#ifdef DEBUG
	fprintf(stderr, "Entering function %s\n", __func__);
#endif
	if(type == EVENT_COAL || type == EVENT_MIGR)
		ev = cache_alloc(cfg->event_cache[type]);
	else
		ev = malloc(evsize[type] + sizeof(int) * 2 * npop_all);

//	ev = (struct event *)GET_OBJ(l);
	ev->type = type;
	ev->t = t;
//	ev->dn = (int *)((char *)l + evsize[type]);
//	ev->sumdn = (int *)((char *)l + evsize[type] + sizeof(int) * npop_all);
	ev->dn_off = evsize[type];
	ev->sumdn_off = evsize[type] + sizeof(int) * npop_all;
	dn_clear(npop_all, GET_DN(ev));

#ifdef DEBUG
	fprintf(stderr, "Allocated event %x at time %.6f with type %d\n", ev, ev->t, ev->type);
#endif
	return ev;
}

void free_event(struct config *cfg, struct event *ev)
{
	if(ev->type == EVENT_COAL || ev->type == EVENT_MIGR)
		cache_free(cfg->event_cache[ev->type], GET_LIST(ev));
	else
		free(GET_LIST(ev));
}

void print_event(struct config *cfg, struct event *ev)
{
	struct list_head *l;
	int i;

	l = GET_LIST(ev);
	fprintf(stderr, "(%x, %x)[type=%d, t=%.6f, dn=(", l, ev, ev->type, ev->t);
	for(i = 0; i < cfg->npop_all; i++)
		fprintf(stderr, "%d, ", GET_DN(ev)[i]);
	fprintf(stderr, ")");
	if(ev->type == EVENT_COAL){
		fprintf(stderr, ", pop=%d", ((struct coal_event *)ev)->pop);
	}else if(ev->type == EVENT_MIGR){
		fprintf(stderr, ", spop=%d, dpop=%d", ((struct migr_event *)ev)->spop, ((struct migr_event *)ev)->dpop);
	}else if(ev->type == EVENT_GROW){
		fprintf(stderr, ", pop=%d, alpha=%.6f", ((struct grow_event *)ev)->pop, ((struct grow_event *)ev)->alpha);
	}else if(ev->type == EVENT_SIZE){	
		fprintf(stderr, ", pop=%d, size=%.6f", ((struct size_event *)ev)->pop, ((struct size_event *)ev)->size);
	}else if(ev->type == EVENT_GMIG){
		fprintf(stderr, ", rmig=%.6f", ((struct gmig_event *)ev)->rmig);
	}else if(ev->type == EVENT_RMIG){
		fprintf(stderr, ", popi=%d, popj=%d, rmig=%.6f", ((struct rmig_event *)ev)->popi, ((struct rmig_event *)ev)->popj, ((struct rmig_event *)ev)->rmig);
	}else if(ev->type == EVENT_GSIZ){
		fprintf(stderr, ", size=%.6f", ((struct gsiz_event *)ev)->size);
	}else if(ev->type == EVENT_GGRO){
		fprintf(stderr, ", alpha=%.6f", ((struct ggro_event *)ev)->alpha);
	}else if(ev->type == EVENT_JOIN){
		fprintf(stderr, ", popi=%d, popj=%d", ((struct join_event *)ev)->popi, ((struct join_event *)ev)->popj);
	}else if(ev->type == EVENT_SPLT){
		fprintf(stderr, ", pop=%d, newpop=%d", ((struct splt_event *)ev)->pop, ((struct splt_event *)ev)->newpop);
	}else if(ev->type == EVENT_SAMP){
	}

	fprintf(stderr, "]");
}

char nucl[] = {'A', 'C', 'G', 'T'};

int nucl_index(int ch)
{
	switch(ch){
		case 'a':case 'A':
			return 0;
		case 'c':case 'C':
			return 1;
		case 'g':case 'G':
			return 2;
		case 't':case 'T':
			return 3;
		default:
			return -1;
	}
}

int fgcompar_r(const void *a, const void *b, int *fgstart)
{
	return fgstart[*(int *)a] - fgstart[*(int *)b];
}

void print_profile(struct profile *prof, FILE *outfp)
{
	int i, j, *fgstart, *fgend, *fgid;
	struct fginfo *fgi;

	fgstart = prof->fgstart;
//	fgend = prof->fgend;
	fgid = prof->fgid;
	fgi = prof->info;

	fprintf(outfp, "%s,%d\n", prof->reffile, prof->chrnum);

	fprintf(outfp, "%d", prof->npop);
	for(i = 0; i < prof->npop; i++)
		fprintf(outfp, ",%d", prof->nfrag_pop[i] + prof->ntrunks[i]);
	fprintf(outfp, "\n");

	for(i = 0; i < prof->nfrag; i++){
		fprintf(outfp, "%d\t%d\t%d\t%d\t%d", fgid[i], fgi[i].pop, fgstart[i], fgi[i].end, fgi[i].nread);
		for(j = 0; j < fgi[i].nread; j++){
			fprintf(outfp, "\t%d\t%d", prof->rdset[i][j].start, prof->rdset[i][j].end);
		}
		fprintf(outfp, "\n");
	}
}

struct profile *generate_profile(char *reffile, int chrnum, int npop, int *nfrags, int fraglen, int paired, int rdlen, int *ntrunks)
{
	struct reference *ref;
	struct profile *prof;
	int len, reflen, nfrag, nread, i, j;
	int *fgstart, *fgend, *fgid, *fgorder;
	struct fginfo *fgi;

	if(paired == 0)
		nread = 1;
	else
		nread = 2;

	prof = malloc(sizeof(struct profile));
	len = strlen(reffile);

	prof->reffile = malloc(sizeof(char) * (len + 1));
	strncpy(prof->reffile, reffile, len);
	prof->reffile[len] = '\0';

	prof->ref = load_reference(reffile);
	prof->chrnum = chrnum;
	prof->npop = npop;

	prof->nfrag_pop = malloc(sizeof(int) * npop);
	memcpy(prof->nfrag_pop, nfrags, sizeof(int) * npop);

	prof->ntrunks = malloc(sizeof(int) * npop);
	memcpy(prof->ntrunks, ntrunks, sizeof(int) * npop);

	nfrag = 0;
	for(i = 0; i < npop; i++)
		nfrag += nfrags[i] + ntrunks[i];

	fgstart = malloc(sizeof(int) * nfrag);
	fgi = malloc(sizeof(struct fginfo) * nfrag);

	reflen = prof->ref->chrlen[chrnum];

	for(i = nfrag = 0; i < npop; i++){
		for(j = 0; j < ntrunks[i]; j++, nfrag++){
			fgi[nfrag].pop = i;
			fgstart[nfrag] = 0;
			fgi[nfrag].end = reflen;
			fgi[nfrag].nread = 1;
			fgi[nfrag].trunk = 0;
		}
	}

	for(i = 0; i < npop; i++){
		for(j = 0; j < nfrags[i]; j++, nfrag++){
			int readpos, fragpos;

			fragpos = dunif(reflen - fraglen);
			readpos = fragpos;
			fgi[nfrag].pop = i;
			fgstart[nfrag] = fragpos;
			fgi[nfrag].end = fragpos + fraglen;
			fgi[nfrag].nread = nread;
			fgi[nfrag].trunk = 0;
		}
	}
	prof->nfrag = nfrag;

	fgorder = malloc(sizeof(int) * nfrag);
	for(i = 0; i < nfrag; i++)
		fgorder[i] = i;

	// Sort fragments by position
	qsort_r(fgorder, nfrag, sizeof(int), fgcompar_r, prof->fgstart);

	prof->fgstart = malloc(sizeof(int) * nfrag);
	prof->fgid = malloc(sizeof(int) * nfrag);
	prof->info = malloc(sizeof(struct fginfo) * nfrag);
	prof->nds = malloc(sizeof(struct sam_node *) * nfrag);

	for(i = 0; i < nfrag; i++){
		int fragpos, readpos, k, nread;

		prof->fgid[i] = i + 1;
		prof->fgstart[i] = fgstart[fgorder[i]];
		prof->info[i] = fgi[fgorder[i]];

		readpos = fragpos = prof->fgstart[i];
		nread = prof->info[i].nread;

		prof->rdset[i] = malloc(nread * sizeof(struct read));
		for(k = 0; k < nread; k++){
			prof->rdset[i][k].start = readpos;
			prof->rdset[i][k].end = readpos + rdlen;
			readpos = fragpos + fraglen - rdlen;
			prof->rdset[i][k].seq = NULL;
			prof->rdset[i][k].qual = NULL;
		}
	}

	free(fgstart);
	free(fgi);
	free(fgorder);

	return prof;
}

/* Destroy object of read profile */
void unload_profile(struct profile *prof)
{
	int i;

	if(prof->rdset){
		for(i = 0; i < prof->nfrag; i++){
			int j;
			for(j = 0; j < prof->info[i].nread; j++){
				if(prof->rdset[i][j].seq)
					free(prof->rdset[i][j].seq);

				if(prof->rdset[i][j].qual)
					free(prof->rdset[i][j].qual);
			}

			free(prof->rdset[i]);
		}
		free(prof->rdset);
	}

	unload_reference(prof->ref);
	free(prof->ntrunks);
	free(prof->nfrag_pop);
	free(prof->fgstart);
	free(prof->fgid);
	free(prof->info);
	free(prof->nds);
	free(prof->reffile);
	free(prof);
}

/* Load read profile. */
struct profile *load_profile(FILE *filp, int genseq)
{
	int npop, nfrag, maxfrag, i, ch, *fgorder, *fgstart, *fgend, *fgid;
	struct profile *prof;
	struct fginfo *fgi;
	struct read **rdset;
	char file[1000];

	/***** Load read profile *****/
	prof = malloc(sizeof(struct profile));
	memset(prof, 0, sizeof(struct profile));

	/* Read reference file name. */
	i = 0;
	while((ch = fgetc(filp)) != ',') file[i++] = ch;
	file[i] = '\0';
	prof->reffile = malloc(sizeof(char) * (i + 1));
	strncpy(prof->reffile, (char *)file, (size_t)(i + 1));
	prof->ref = load_reference(prof->reffile);
	ch = read_integer(filp, &prof->chrnum);

	/* Read number of subpopulations. */
	ch = read_integer(filp, &npop);
	prof->npop = npop;
	if(ch != ','){
		fprintf(stderr, "Expect for , after number of subpopulations.\n");
		goto abnormal;
	}

	prof->nfrag_pop = malloc(sizeof(int) * npop);
	nfrag = 0;
	for(i = 0; i < npop; i++){
		/* Read number of fragments. */
		ch = read_integer(filp, &prof->nfrag_pop[i]);
		nfrag += prof->nfrag_pop[i];
		if(ch != ',' && ch != '\n'){
			fprintf(stderr, "Expect ',' or newline in line 2.\n");
			goto abnormal;
		}
	}
	if(i < npop){
		fprintf(stderr, "Too few number of subpopulations.\n");
		goto abnormal;
	}
	prof->nfrag = nfrag;
	while(ch != '\n') ch = fgetc(filp);

	if(ch != '\n'){
		fprintf(stderr, "Expect for , after read numbers.\n");
		goto abnormal;
	}

	/* Load fragment information from input file. */
	fgstart = malloc(sizeof(int) * (nfrag + 1));
	fgid = malloc(sizeof(int) * (nfrag + 1));
	fgi = malloc(sizeof(struct fginfo) * (nfrag + 1));

	if(genseq)
		rdset = malloc(sizeof(struct read *) * (nfrag + 1));
	fgorder = malloc(sizeof(int) * (nfrag + 1));

	for(i = 0; i < nfrag; i++){
		struct frag *fg;
		int j, pop, nread, seqlen, intbuf;

		fgorder[i] = i;

		ch = read_integer(filp, &fgid[i]);
		ch = read_integer(filp, &pop);
		ch = read_integer(filp, &fgstart[i]);
		ch = read_integer(filp, &fgi[i].end);
		ch = read_integer(filp, &nread);
		fgi[i].pop = pop;
		fgi[i].nread = nread;
		fgi[i].trunk = 0;

		if(genseq){
			rdset[i] = malloc(sizeof(struct read) * fgi[i].nread);

			for(j = 0; j < fgi[i].nread; j++){
				ch = read_integer(filp, &rdset[i][j].start);
				ch = read_integer(filp, &rdset[i][j].end);

				seqlen = rdset[i][j].end - rdset[i][j].start;
				rdset[i][j].seq = malloc(sizeof(char) * (seqlen + 1));
				memset(rdset[i][j].seq, 0, sizeof(char) * (seqlen + 1));
				rdset[i][j].qual = NULL;
			}
		}else{
			for(j = 0; j < fgi[i].nread; j++){
				ch = read_integer(filp, &intbuf);
				ch = read_integer(filp, &intbuf);
			}
		}

		while(ch != '\n') ch = fgetc(filp);
	}

	fgstart[nfrag] = fgi[nfrag].end = INT_MAX;
	qsort_r((void *)fgorder, (size_t)nfrag, sizeof(int), fgcompar_r, fgstart);
	prof->fgstart = malloc(sizeof(int) * (nfrag + 1));
	prof->fgid = malloc(sizeof(int) * (nfrag + 1));
	prof->info = malloc(sizeof(struct fginfo) * (nfrag + 1));

	for(i = 0; i < nfrag; i++){
		prof->fgstart[i] = fgstart[fgorder[i]];
		prof->fgid[i] = fgid[fgorder[i]];
		prof->info[i] = fgi[fgorder[i]];
	}

	if(genseq){
		prof->rdset = malloc(sizeof(struct read *) * (nfrag + 1));
		for(i = 0; i < nfrag; i++)
			prof->rdset[i] = rdset[fgorder[i]];
	}

	prof->nds = malloc(sizeof(struct sam_node *) * (nfrag + 1));

	free(fgstart);
	free(fgid);
	free(fgi);

	if(genseq)
		free(rdset);
	free(fgorder);

finish:
	return prof;

abnormal:
	if(prof->ref)
		unload_reference(prof->ref);

	if(prof->reffile)
		free(prof->reffile);

	if(prof->nfrag_pop)
		free(prof->nfrag_pop);

	free(prof);
	prof = NULL;
	goto finish;
}

void print_fragment(FILE *outfp, int fgstart, int fgid, struct fginfo *fgi, struct read *rdset)
{
	struct read *rd;
	int i;
#ifdef DEBUG
	fprintf(stderr, "Entering %s: fgid=%x, nread=%d\n", __func__, fgid, fgi->nread);
#endif
	for(i = 0; i < fgi->nread; i++){
		rd = &rdset[i];
		fprintf(outfp, ">%d:%d %d %d-%d %d-%d\n", fgid, i + 1, fgi->pop, fgstart, fgi->end, rd->start, rd->end);
		fprintf(outfp, "%s\n", rd->seq);
	}
#ifdef DEBUG
	fprintf(stderr, "Finishing %s\n", __func__);
#endif
}

/* seed: Seed of random number generator.
   print_tree: This is set either 0 or 1. If the value is 0, the program output reads and do not print local genealogies. If the value is 1, local genealogies are printed and read sequence are not printed.
   maxfrag: Maximum number of reads in a group.
   theta: Mutation rates of each populaion.
   mmig: Migration rate matrix between populations that exists at time 0. If there is npop populations, this array should be of length npop * npop. Give NULL if no migration are allowed.
   grate: Growth rates of populations that exists at time 0. Give NULL if population size does not change over time.
   */
struct config *create_config(int seed, int print_tree, int gensam, FILE *treefp, FILE *readfp, int maxfrag, struct profile *prof, double rho)
{
	struct config *cfg;
	struct event *ev;
	int npop, nsplt;
	double *mmig;
	int i, j, reflen;

	cfg = malloc(sizeof(struct config));
	memset(cfg, 0, sizeof(struct config));

	cfg->prof = prof;

	/* Set up basic configuration */
	cfg->npop = cfg->npop_all = npop = prof->npop;
	cfg->nsplt = nsplt = 0;

	cfg->seed = seed;
	cfg->print_tree = print_tree;
	cfg->gensam = gensam;
	cfg->treefp = treefp;
	cfg->readfp = readfp;
	cfg->maxfrag = maxfrag;
	cfg->rho = rho;

	cfg->mmut = malloc(sizeof(struct mutation) * cfg->npop);
	memset(cfg->mmut, 0, sizeof(struct mutation) * cfg->npop);

	cfg->size = malloc(prof->npop * sizeof(double));
	for(i = 0; i < prof->npop; i++)
		cfg->size[i] = 1;

	cfg->grate = malloc(sizeof(double) * npop);
	memset(cfg->grate, 0, sizeof(double) * npop);

	/* Set up migration matrix at time 0 */
	cfg->mmig = malloc(sizeof(double *) * cfg->npop);
	mmig = malloc(sizeof(double) * cfg->npop * cfg->npop);
	for(i = 0; i < cfg->npop; i++){
		cfg->mmig[i] = mmig;
		for(j = 0; j < cfg->npop; j++){
			if(i == j)
				mmig[j] = 1;
			else
				mmig[j] = (double)1 / (cfg->npop - 1);
		}
		mmig += cfg->npop;
	}

	if(cfg->npop == 1){
		cfg->mmig[0][0] = 0;
	}

	/* Set up empty event list (contains only one dummy event) */
	ev = alloc_event(cfg, EVENT_GSIZ, INFINITY);

	list_init(&cfg->evlist);
	cfg->ndevents = 0;
	cfg->devents = malloc(sizeof(struct event *) * 10);

	list_append(&cfg->evlist, ev);
	((struct gsiz_event *)ev)->size = 1;

	reflen = prof->ref->chrlen[prof->chrnum];

	return cfg;
}

void dump_config(struct config *cfg)
{
	struct list_head *l;
	int i, j;

	fprintf(stderr, "Starting %s:\n", __func__);
	fprintf(stderr, " seed=%d, print_tree=%d, treefp=%x, rho=%.6f, npop=%d, nsplt=%d, npop_all=%d, maxfrag=%d\n",
			cfg->seed, cfg->print_tree, cfg->treefp, cfg->rho, cfg->npop, cfg->nsplt, cfg->npop_all, cfg->maxfrag);

	fprintf(stderr, " size: [");
	for(i = 0; i < cfg->npop; i++){
		fprintf(stderr, "%.6f, ", cfg->size[i]);
	}
	fprintf(stderr, "]\n");

	fprintf(stderr, " mmig: [\n");
	for(i = 0; i < cfg->npop; i++){
		fprintf(stderr, "       [");
		for(j = 0; j < cfg->npop; j++){
			fprintf(stderr, "%.6f, ", cfg->mmig[i][j]);
		}

		fprintf(stderr, "]\n");
	}
	fprintf(stderr, " ]\n");

	fprintf(stderr, " grate: [");
	for(i = 0; i < cfg->npop; i++){
		fprintf(stderr, "%.6f, ", cfg->grate[i]);
	}
	fprintf(stderr, "]\n");

	// Print mutation models
	for(i = 0; i < cfg->npop_all; i++){
		if(cfg->mmut[i])
			dump_mutation_model(cfg->mmut[i]);
		else
			fprintf(stderr, "Mutation model of population %d is not yet registered.\n", i + 1);
	}

	// Print event list
	fprintf(stderr, "\nPrint event list:\n");
	l = cfg->evlist.front;
	while(l){
		print_event(cfg, (struct event *)GET_OBJ(l));
		l = l->next;
		fprintf(stderr, "\n");
	}
	fprintf(stderr, "\n");
}

int check_config(struct config *cfg)
{
	return 0;
}

void destroy_config(struct config *cfg)
{
	struct list_head *l, *tmp;
	struct event *ev;
	int i;

	free(cfg->evlist.front);
	list_init(&cfg->evlist);

	for(i = 0; i < cfg->ndevents; i++){
		if(cfg->devents[i]->type == EVENT_SPLT)
			node_set_destroy(&((struct splt_event *)cfg->devents[i])->ndls);
		else if(cfg->devents[i]->type == EVENT_JOIN)
			node_set_destroy(&((struct join_event *)cfg->devents[i])->ndls);
		l = GET_LIST(cfg->devents[i]);
		free(l);
	}

	free(cfg->devents);
	free(cfg->mmut);
	free(cfg->grate);
	free(cfg->mmig[0]);
	free(cfg->mmig);
	free(cfg->size);

	free(cfg);
}

int register_mutation_model(struct config *cfg, int pop, struct mutation *mmut)
{
	int i;

	cfg->mmut[pop] = mmut;

	return 0;

}

int set_growth_rates(struct config *cfg, double *grate)
{
	int i;

	for(i = 0; i < cfg->npop; i++)
		cfg->grate[i] = grate[i];

	return 0;
}

int set_migration_matrix(struct config *cfg, double *mmig)
{
	double *ptr;
	int j, k;

	ptr = mmig;
	for(j = 0; j < cfg->npop; j++){
		memset(cfg->mmig[j], 0, cfg->npop * sizeof(double));
		for(k = 0; k < cfg->npop; k++){
			if(j != k){
				cfg->mmig[j][k] = mmig[j * cfg->npop + k];
				cfg->mmig[j][j] += mmig[j * cfg->npop + k];
			}
		}
	}


	return 0;
}

int add_event_ggro(struct config *cfg, double t, double alpha)
{
	struct list_head *evl;
	struct event *ev;

	ev = alloc_event(cfg, EVENT_GGRO, t);
	((struct ggro_event *)ev)->alpha = alpha;

	evl = cfg->evlist.front;
	while(evl && ((struct event *)GET_OBJ(evl))->t < t) evl = evl->next;
	list_insbefore(evl, ev);
	cfg->ndevents++;

	return 0;
}

int add_event_grow(struct config *cfg, double t, int pop, double alpha)
{
	struct list_head *evl;
	struct event *ev;

	ev = alloc_event(cfg, EVENT_GROW, t);

	((struct grow_event *)ev)->pop = pop;
	((struct grow_event *)ev)->alpha = alpha;

	evl = cfg->evlist.front;
	while(evl && ((struct event *)GET_OBJ(evl))->t < t) evl = evl->next;
	list_insbefore(evl, ev);
	cfg->ndevents++;

	return 0;
}

int add_event_gmig(struct config *cfg, double t, double rmig)
{
	struct list_head *evl;
	struct event *ev;

	ev = alloc_event(cfg, EVENT_GMIG, t);
	((struct gmig_event *)ev)->rmig = rmig;

	evl = cfg->evlist.front;
	while(evl && ((struct event *)GET_OBJ(evl))->t < t) evl = evl->next;
	list_insbefore(evl, ev);
	cfg->ndevents++;

	return 0;
}

int add_event_rmig(struct config *cfg, double t, int popi, int popj, double rmig)
{
	struct list_head *evl;
	struct event *ev;

	ev = alloc_event(cfg, EVENT_RMIG, t);

	((struct rmig_event *)ev)->popi = popi;
	((struct rmig_event *)ev)->popj = popj;
	((struct rmig_event *)ev)->rmig = rmig;

	evl = cfg->evlist.front;
	while(evl && ((struct event *)GET_OBJ(evl))->t < t) evl = evl->next;
	list_insbefore(evl, ev);
	cfg->ndevents++;

	return 0;
}

int add_event_gsiz(struct config *cfg, double t, double size)
{
	struct list_head *evl;
	struct event *ev;

	ev = alloc_event(cfg, EVENT_GSIZ, t);
	((struct gsiz_event *)ev)->size = size;

	evl = cfg->evlist.front;
	while(evl && ((struct event *)GET_OBJ(evl))->t < t) evl = evl->next;
	list_insbefore(evl, ev);
	cfg->ndevents++;

	return 0;
}

int add_event_join(struct config *cfg, double t, int popi, int popj)
{
	struct list_head *evl;
	struct event *ev;

	ev = alloc_event(cfg, EVENT_JOIN, t);
	((struct join_event *)ev)->popi = popi;
	((struct join_event *)ev)->popj = popj;
	node_set_init(&((struct join_event *)ev)->ndls, cfg->maxfrag);

	evl = cfg->evlist.front;
	while(evl && ((struct event *)GET_OBJ(evl))->t < t) evl = evl->next;
	list_insbefore(evl, ev);
	cfg->ndevents++;

	return 0;
}

int add_event_splt(struct config *cfg, double t, int pop, double prop)
{
	struct list_head *evl;
	struct event *ev;
	int npop_all;

	npop_all = cfg->npop + cfg->nsplt;
	ev = alloc_event(cfg, EVENT_SPLT, t);
	((struct splt_event *)ev)->pop = pop;
	((struct splt_event *)ev)->newpop = cfg->npop_all++;
	((struct splt_event *)ev)->prop = prop;
	node_set_init(&((struct splt_event *)ev)->ndls, cfg->maxfrag);

	evl = cfg->evlist.front;
	while(evl && ((struct event *)GET_OBJ(evl))->t < t) evl = evl->next;
	list_insbefore(evl, ev);
	cfg->mmut = realloc(cfg->mmut, sizeof(struct mutation *) * npop_all);
	cfg->mmut[cfg->npop_all - 1] = NULL;
	cfg->size = realloc(cfg->size, sizeof(double) * npop_all);
	cfg->size[cfg->npop_all - 1] = cfg->size[pop];
	cfg->grate = realloc(cfg->grate, sizeof(double) * npop_all);
	cfg->grate[cfg->npop_all - 1] = 0;
	cfg->ndevents++;

	return 0;
}

int add_event_size(struct config *cfg, double t, int pop, double size)
{
	struct list_head *evl;
	struct event *ev;

	ev = alloc_event(cfg, EVENT_SIZE, t);
	((struct size_event *)ev)->pop = pop;
	((struct size_event *)ev)->size = size;

	evl = cfg->evlist.front;
	while(evl && ((struct event *)GET_OBJ(evl))->t < t) evl = evl->next;
	list_insbefore(evl, ev);
	cfg->ndevents++;

	return 0;
}

int add_event_samp(struct config *cfg, double t, int pop, double size)
{
	struct list_head *evl;
	struct event *ev;

	ev = alloc_event(cfg, EVENT_SAMP, t);

	evl = cfg->evlist.front;
	while(evl && ((struct event *)GET_OBJ(evl))->t < t) evl = evl->next;
	list_insbefore(evl, ev);
	cfg->ndevents++;

	return 0;
}

struct reference *load_reference(char *file)
{
	struct reference *ref;
	int ch, nchr, maxchr;

	maxchr = 40;
	nchr = 0;

	ref = malloc(sizeof(struct reference));
	ref->seq_start = malloc(sizeof(int) * maxchr);
	ref->chrlen = malloc(sizeof(int) * maxchr);

	/* Skip header of reference genome. */
	ref->filp = fopen(file, "r");
	ch = fgetc(ref->filp);
	while(ch == '>'){
		double chrlen;

		// Load chromosome header
		while((ch = fgetc(ref->filp)) != EOF && ch != '\n');
		if(ch == EOF){
			fprintf(stderr, "No sequence in reference file.");
		}

		// Load chromosome sequence
		ref->seq_start[nchr] = ftell(ref->filp);
		chrlen = 0;
		ch = fgetc(ref->filp);
		while(ch != EOF && ch != '>'){
			chrlen++;
			while((ch = fgetc(ref->filp)) != '\n') chrlen++;
			ch = fgetc(ref->filp);
		}
		ref->chrlen[nchr++] = chrlen;;
	}
	ref->nchr = nchr;

	return ref;
}

void unload_reference(struct reference *ref)
{
	free(ref->seq_start);
	free(ref->chrlen);
	fclose(ref->filp);
	free(ref);
}

void reload_reference(struct reference *ref, int chrnum)
{
	fseek(ref->filp, ref->seq_start[chrnum], SEEK_SET);
	ref->curr = 0;
}

void load_chr(struct reference *ref, int chrnum, char **strp)
{
	int reflen, i, curr, ch;
	char *str;

	reload_reference(ref, chrnum);
	reflen = ref->chrlen[chrnum];

	str = malloc(sizeof(char) * (reflen + 1));
	curr = 0;
	while(curr < reflen){
		while(!isalpha(ch = fgetc(ref->filp)) && ch > 0) curr++;
		str[curr++] = ch;
	}

	str[curr] = '\0';
	*strp = str;
}

