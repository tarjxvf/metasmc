#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include "util.h"
#include "global.h"
#include "smc.h"
#include "mutation.h"
#include "evindex.h"

void list_init(struct list *ls)
{
	ls->front = NULL;
	ls->rear = &ls->front;
	ls->n = 0;
}

void list_concat(struct list *dst, struct list *src)
{
	*dst->rear = src->front;
	src->front->prev = dst->rear;
	dst->rear = src->rear;
}

/* Add an item in the front of the list. */
void list_add(struct list *ls, void *item)
{
	struct list_head *l;

	l = GET_LIST(item);
	__list_add(ls, l);
}

/* Insert an item before an item. */
void list_insbefore(struct list_head *ref, void *item)
{
	struct list_head *l;

	l = GET_LIST(item);
	__list_insbefore(ref, l);
}

/* Append an item after a  list */
void list_append(struct list *ls, void *item)
{
	struct list_head *l;

	l = GET_LIST(item);
	__list_append(ls, l);
}

/* Remove an item from the list. */
void list_remove(struct list *ls, void *item)
{
	struct list_head *l;
	l = GET_LIST(item);
	__list_remove(ls, l);
}

void list_print(struct list_head **head)
{
	struct list_head **ptr;
	ptr = head;

	printf("head=%x", ptr);
	while(*ptr){
		printf("->[ptr=%x, *ptr=%x, prev=%x, next=%x, &next=%x]", ptr, *ptr, (*ptr)->prev, (*ptr)->next, &(*ptr)->next);
		ptr = &(*ptr)->next;
	}
	printf("\n");
}

/* Add an item in the front of the list. *
void list_add(struct list_head **head, struct list_head *item)
{
	if(*head)
		(*head)->prev = &item->next;
	item->next = *head;
	*head = item;
	item->prev = head;
}

/* Insert an item before an item. *
void list_insbefore(struct list_head *ref, struct list_head *item)
{
	item->next = ref;
	item->prev = ref->prev;
	ref->prev = &item->next;
	*item->prev = item;
}

/* Append an item after a  list *
void list_append(struct list_head **head, struct list_head *item)
{
	struct list_head **ptr;

	ptr = head;
	while(*ptr)
		ptr = &(*ptr)->next;
	*ptr = item;
//	item->next = NULL;
	item->prev = ptr;
}

/* Remove an item from the list. *
void list_remove(struct list_head *item)
{
	*item->prev = item->next;
	if(item->next)
		item->next->prev = item->prev;
}

void list_print(struct list_head **head)
{
	struct list_head **ptr;
	ptr = head;

	printf("head=%x", ptr);
	while(*ptr){
		printf("->[ptr=%x, *ptr=%x, prev=%x, next=%x, &next=%x]", ptr, *ptr, (*ptr)->prev, (*ptr)->next, &(*ptr)->next);
		ptr = &(*ptr)->next;
	}
	printf("\n");
}*/

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

void add_event(struct genealogy *G, struct event *ev)
{
	evindex_insert(G->evidx, ev);
}

void remove_event(struct genealogy *G, struct event *ev)
{
	struct list_head *l;

#ifdef DEBUG
	fprintf(stderr, "Entering function %s, event=%x, type=%d, t=%.6f\n", __func__, ev, ev->type, ev->t);
#endif
	if(ev->type == EVENT_JOIN || ev->type == EVENT_SPLT){
		struct rb_traverser tr;

//		rb_t_find(&tr, G->evidx->idx->tree, ev);

		if(ev->type == EVENT_JOIN){
			struct join_event *jev;

			jev = (struct join_event *)ev;
			memset(G->evidx->dn, 0, sizeof(int) * G->cfg->npop_all);
			G->evidx->dn[jev->popi] = -1;
			G->evidx->dn[jev->popj] = 1;

			jev->dn[jev->popi]--;
			jev->dn[jev->popj]++;

		}else if(ev->type == EVENT_SPLT){
			struct splt_event *sev;

			sev = (struct splt_event *)ev;
			memset(G->evidx->dn, 0, sizeof(int) * G->cfg->npop_all);
			G->evidx->dn[sev->newpop] = -1;
			G->evidx->dn[sev->pop] = 1;

			sev->dn[sev->pop]++;
			sev->dn[sev->newpop]--;
		}

		// Update tree values
//		evindex_propgate(&tr, G->cfg->npop_all, G->evidx->dn);

	}else{
		evindex_delete(G->evidx, ev);
		l = GET_LIST(ev);
		free(l);
	}
}

size_t evsize[] = {sizeof(struct list_head) + sizeof(struct coal_event), sizeof(struct list_head) + sizeof(struct migr_event), sizeof(struct list_head) + sizeof(struct grow_event), sizeof(struct list_head) + sizeof(struct size_event), sizeof(struct list_head) + sizeof(struct rmig_event), sizeof(struct list_head) + sizeof(struct gmig_event), sizeof(struct list_head) + sizeof(struct gsiz_event), sizeof(struct list_head) + sizeof(struct ggro_event), sizeof(struct list_head) + sizeof(struct join_event), sizeof(struct list_head) + sizeof(struct splt_event), sizeof(struct list_head) + sizeof(struct event), sizeof(struct list_head) + sizeof(struct event), sizeof(struct list_head) + sizeof(struct samp_event)};

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
	l = malloc(evsize[type] + sizeof(int) * 2 * npop_all);
	ev = (struct event *)GET_OBJ(l);
	ev->type = type;
	ev->t = t;
	ev->dn = (int *)((char *)l + evsize[type]);
	ev->sumdn = (int *)((char *)l + evsize[type] + sizeof(int) * npop_all);
	memset(ev->dn, 0, sizeof(int) * 2 * npop_all);

#ifdef DEBUG
	fprintf(stderr, "Allocated event %x at time %.6f with type %d\n", ev, ev->t, ev->type);
#endif
	return ev;
}

void print_event(struct config *cfg, struct event *ev)
{
	struct list_head *l;
	int i;

	l = GET_LIST(ev);
	fprintf(stderr, "(%x, %x)[type=%d, t=%.6f, dn=(", l, ev, ev->type, ev->t);
	for(i = 0; i < cfg->npop_all; i++)
		fprintf(stderr, "%d, ", ev->dn[i]);
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
/*	}else if(ev->type == EVENT_MMUT){
		fprintf(stderr, ", pop=%d", ((struct mmut_event *)ev)->pop);
		dump_mutation_model(((struct mmut_event *)ev)->mmut);*/
	}else if(ev->type == EVENT_SAMP){
	}

	fprintf(stderr, "]");
}

int fgcompar(const void *a, const void *b)
{
	return ((struct frag *)a)->start - ((struct frag *)b)->start;
}

void print_profile(struct profile *prof, FILE *outfp)
{
	int i, j;

	fprintf(outfp, "%s,%d\n", prof->reffile, prof->chrnum);

	fprintf(outfp, "%d", prof->npop);
	for(i = 0; i < prof->npop; i++)
		fprintf(outfp, ",%d", prof->nfrag_pop[i] + prof->ntrunks[i]);
	fprintf(outfp, "\n");

	for(i = 0; i < prof->nfrag; i++){
		fprintf(outfp, "%d\t%d\t%d\t%d\t%d", prof->fgset[i].id, prof->fgset[i].pop, prof->fgset[i].start, prof->fgset[i].end, prof->fgset[i].nread);
		for(j = 0; j < prof->fgset[i].nread; j++){
			fprintf(outfp, "\t%d\t%d", prof->fgset[i].rd[j].start, prof->fgset[i].rd[j].end);
		}
		fprintf(outfp, "\n");
	}
}

struct profile *generate_profile(char *reffile, int chrnum, int npop, int *nfrags, int fraglen, int paired, int rdlen, int *ntrunks)
{
	struct reference *ref;
	struct profile *prof;
	int len, reflen, nfrag, nread, i, j;

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
	prof->fgset = malloc(sizeof(struct frag) * nfrag);

	reflen = prof->ref->chrlen[chrnum];

	for(i = nfrag = 0; i < npop; i++){
		for(j = 0; j < ntrunks[i]; j++, nfrag++){
			prof->fgset[nfrag].pop = i;
			prof->fgset[nfrag].start = 0;
			prof->fgset[nfrag].end = reflen;
			prof->fgset[nfrag].nread = 1;
			prof->fgset[nfrag].rd = malloc(sizeof(struct read) * 1);
			prof->fgset[nfrag].trunk = 0;

			prof->fgset[nfrag].rd[0].start = 0;
			prof->fgset[nfrag].rd[0].end = reflen;
			prof->fgset[nfrag].rd[0].seq = NULL;
			prof->fgset[nfrag].rd[0].qual = NULL;
		}
	}

	for(i = 0; i < npop; i++){
		for(j = 0; j < nfrags[i]; j++, nfrag++){
			int fragpos, readpos, k;

			fragpos = dunif(reflen - fraglen);
			readpos = fragpos;
			prof->fgset[nfrag].pop = i;
			prof->fgset[nfrag].start = fragpos;
			prof->fgset[nfrag].end = fragpos + fraglen;
			prof->fgset[nfrag].nread = nread;
			prof->fgset[nfrag].rd = malloc(sizeof(struct read) * nread);
			prof->fgset[nfrag].trunk = 0;

			readpos = fragpos;
			for(k = 0; k < nread; k++){
				prof->fgset[nfrag].rd[k].start = readpos;
				prof->fgset[nfrag].rd[k].end = readpos + rdlen;
				readpos = fragpos + fraglen - rdlen;
				prof->fgset[nfrag].rd[k].seq = NULL;
				prof->fgset[nfrag].rd[k].qual = NULL;
			}
		}
	}
	prof->nfrag = nfrag;

	// Sort fragments by position
	qsort(prof->fgset, nfrag, sizeof(struct frag), fgcompar);
	for(i = 0; i < nfrag; i++)
		prof->fgset[i].id = i + 1;
//	unload_reference(prof->ref);

	return prof;
}

/* Destroy object of read profile */
void unload_profile(struct profile *prof)
{
	int i;

	for(i = 0; i < prof->nfrag; i++){
		int j;
		for(j = 0; j < prof->fgset[i].nread; j++){
			if(prof->fgset[i].rd[j].seq)
				free(prof->fgset[i].rd[j].seq);

			if(prof->fgset[i].rd[j].qual)
				free(prof->fgset[i].rd[j].qual);
		}

		free(prof->fgset[i].rd);
	}
	unload_reference(prof->ref);
	free(prof->ntrunks);
	free(prof->nfrag_pop);
//	free(prof->reflen);
	free(prof->fgset);
	free(prof->reffile);
//	free(prof->popsize);
	free(prof);
}

/* Load read profile. */
struct profile *load_profile(FILE *filp)
{
	int npop, nfrag, maxfrag, i, ch;
	struct profile *prof;
	struct frag *fgset;
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

	/* Read length of reference genome. */
/*	ch = read_integer(filp, &prof->reflen);
	if(ch != ','){
		fprintf(stderr, "Expect for , after reference genome length.\n");
		goto abnormal;
	}*/

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

	/* Read number of subpopulations. */
//	ch = read_integer(filp, &npop);
//	prof->npop = npop;

	/* Read default population size. */
//	prof->popsize = malloc(sizeof(double) * npop);
//	prof->popsize[0] = 1;
//	for(i = 1; i < npop; i++)
//	for(i = 0; i < npop; i++)
//		ch = read_double(filp, &prof->popsize[i]);

	if(ch != '\n'){
		fprintf(stderr, "Expect for , after read numbers.\n");
		goto abnormal;
	}

	/* Load fragment information from input file. */
	fgset = malloc(sizeof(struct frag) * (nfrag + 1));
	memset(fgset, 0, sizeof(struct frag) * (nfrag + 1));
	for(i = 0; i < nfrag; i++){
		struct frag *fg;
		int j;

		fg = &fgset[i];

		ch = read_integer(filp, &fg->id);
		ch = read_integer(filp, &fg->pop);
//		ch = read_integer(filp, &fg->chr);
//		ch = read_double(filp, &fg->start);
//		ch = read_double(filp, &fg->end);
		ch = read_integer(filp, &fg->start);
		ch = read_integer(filp, &fg->end);
		ch = read_integer(filp, &fg->nread);

		fg->rd = malloc(sizeof(struct read) * fg->nread);
		memset(fg->rd, 0, sizeof(struct read) * fg->nread);
		fg->trunk = 0;

		for(j = 0; j < fg->nread; j++){
			int seqlen;

			ch = read_integer(filp, &fg->rd[j].start);
			ch = read_integer(filp, &fg->rd[j].end);

			seqlen = fg->rd[j].end - fg->rd[j].start;
			fg->rd[j].seq = malloc(sizeof(char) * (seqlen + 1));
			memset(fg->rd[j].seq, 0, sizeof(char) * (seqlen + 1));
		}

		while(ch != '\n') ch = fgetc(filp);
	}

	fgset[nfrag].start = fgset[nfrag].end = INT_MAX;
	prof->fgset = fgset;

finish:
	return prof;

abnormal:
	if(prof->ref)
		unload_reference(prof->ref);

	if(prof->reffile)
		free(prof->reffile);

	if(prof->nfrag_pop)
		free(prof->nfrag_pop);
//	if(prof->popsize)
//		free(prof->popsize);

	free(prof);
	prof = NULL;
	goto finish;
}

void print_fragment(FILE *outfp, struct frag *fg)
{
	int i;
#ifdef DEBUG
	fprintf(stderr, "Entering %s: fg=%x, nread=%d\n", __func__, fg, fg->nread);
#endif
	for(i = 0; i < fg->nread; i++){
		struct read *rd;

		rd = &fg->rd[i];
		fprintf(outfp, ">%d:%d %d %d-%d %d-%d\n", fg->id, i + 1, fg->pop, fg->start, fg->end, rd->start, rd->end);
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
	int i, j;

	cfg = malloc(sizeof(struct config));
	memset(cfg, 0, sizeof(struct config));

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
//	memcpy(cfg->size, prof->popsize, prof->npop * sizeof(double));

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

	/* Initialize object caches. */
/*	cfg->node_cache[NODE_COAL] = mem_cache_create(sizeof(struct coal_node), NULL, NULL);
	cfg->node_cache[NODE_MIGR] = mem_cache_create(sizeof(struct migr_node), NULL, NULL);
	cfg->node_cache[NODE_XOVER] = mem_cache_create(sizeof(struct xover_node), NULL, NULL);
	cfg->node_cache[NODE_SAM] = mem_cache_create(sizeof(struct list_head) + sizeof(struct sam_node), NULL, NULL);
	cfg->node_cache[NODE_FLOAT] = mem_cache_create(sizeof(struct dummy_node), NULL, NULL);

	cfg->event_cache[EVENT_COAL] = mem_cache_create(sizeof(struct list_head) + sizeof(struct coal_event), NULL, NULL);
	cfg->event_cache[EVENT_MIGR] = mem_cache_create(sizeof(struct list_head) + sizeof(struct migr_event), NULL, NULL);

	cfg->edge_cache = mem_cache_create(sizeof(struct list_head) + sizeof(struct edge), NULL, NULL);*/
	cfg->frag_cache = mem_cache_create(sizeof(struct list_head) + sizeof(struct frag *), NULL, NULL);

	/* Set up empty event list (contains only one dummy event) */
	ev = alloc_event(cfg, EVENT_GSIZ, INFINITY);

	list_init(&cfg->evlist);
	cfg->ndevents = 0;
	cfg->devents = malloc(sizeof(struct event *) * 10);

	list_append(&cfg->evlist, ev);
	((struct gsiz_event *)ev)->size = 1;

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

	// Release space of event list
//	l = cfg->evlist.front;
/*	while(l){
		tmp = l;
		ev = (struct event *)GET_OBJ(l);
		l = l->next;
//		if(ev->type == EVENT_MMUT)
//			free(((struct mmut_event *)ev)->mmut);
		free(tmp);
	}*/
	free(cfg->evlist.front);
	list_init(&cfg->evlist);
	for(i = 0; i < cfg->ndevents; i++){
		l = GET_LIST(cfg->devents[i]);
		free(l);
	}

/*	mem_cache_destroy(cfg->node_cache[NODE_COAL]);
	mem_cache_destroy(cfg->node_cache[NODE_MIGR]);
	mem_cache_destroy(cfg->node_cache[NODE_XOVER]);
	mem_cache_destroy(cfg->node_cache[NODE_SAM]);
	mem_cache_destroy(cfg->node_cache[NODE_FLOAT]);

	mem_cache_destroy(cfg->event_cache[EVENT_COAL]);
	mem_cache_destroy(cfg->event_cache[EVENT_MIGR]);

	mem_cache_destroy(cfg->edge_cache);*/
	if(cfg->frag_cache)
		mem_cache_destroy(cfg->frag_cache);

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

//	for(i = 0; i < cfg->npop; i++)
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
///	cfg->nsplt++;
	((struct splt_event *)ev)->prop = prop;

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

/*int add_event_mmut(struct config *cfg, double t, int pop, double theta, int model, double *pars, double *pi)
{
	struct list_head *evl;
	struct mmut_event *mev;
	struct event *ev;

	ev = alloc_event(cfg, EVENT_MMUT, t);
	mev = (struct mmut_event *)ev;
	mev->pop = pop;

	evl = cfg->evlist.front;
	while(evl && ((struct event *)GET_OBJ(evl))->t < t) evl = evl->next;
	list_insbefore(evl, ev);

	mev->mmut = malloc(sizeof(struct mutation));
	mev->mmut->theta = theta;
	mev->mmut->model = model;
	memcpy(mev->mmut->mpar, pars, sizeof(double) * npar[model]);
	memcpy(mev->mmut->pi, pi, sizeof(double) * NUM_NUCS);
	init_mutation_model(((struct mmut_event *)ev)->mmut);

	return 0;
}*/

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
//		ref->curr = 0;
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

