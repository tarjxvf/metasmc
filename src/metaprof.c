/* This program load a reference genome database and a abundance profile, then produce read profiles of each species in the abundance profile */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "rand.h"
#include "util.h"
#include "global.h"

struct dbentry {
	int taxid;
	char *name;
	char *reffile;
	int copy;
};

struct database {
	int nline;
	struct dbentry *dbe;
};

struct abund {
	int tax;
	double abund;
	int nfrag;
};

struct abdfile {
	int npop;
	struct abund *abd;
};

void usage(char *prog)
{
	fprintf(stderr, "%s -d dbpath -a abpath -o outdir -x prefix -P <platform> -f fraglen -n <Number of fragments> [-p -r readlen -t ntrunk]\n", prog);
	fprintf(stderr, "-d dbpath: Path to database index file.\n");
	fprintf(stderr, "-a abpath: Path to abundance profile.\n");
	fprintf(stderr, "-o outdir: Output directory.\n");
	fprintf(stderr, "-x prefix: Prefix of output file.\n");
	fprintf(stderr, "-f fraglen: Length of the fragments.\n");
	fprintf(stderr, "-n nfrag: Total numbr of fragments. The number of fragment of each subpopulation is drawn from multinomial distribution.\n");
	fprintf(stderr, "-p: Generate paired-end reads.\n");
	fprintf(stderr, "-r readlen: If -p is open, this option must be used to specify length of paired-end reads.\n");
	fprintf(stderr, "-t ntrunk: Number of complete chromosomes in trunk genealogy.\n");
}

int dbcompar(struct dbentry *a, struct dbentry *b)
{
	return a->taxid - b->taxid;
}

int read_string(FILE *fp, char **strp, char sep)
{
	int maxch, i, ch;
	char *str;

	i = maxch = 0;
	str = NULL;
	do{
		maxch += 100;
		str = realloc(str, sizeof(char) * maxch);
		while(i < maxch && (ch = fgetc(fp)) != sep && ch != '\n' && ch > 0) str[i++] = ch;
	}while(i >= maxch);
	str[i] = '\0';
	str = realloc(str, sizeof(char) * (i + 1));
	*strp = str;

	return ch;
}

struct database *load_database(char *dbfile)
{
	struct dbentry *dbe;
	struct database *db;
	int maxentry, nline, ch, i;
	FILE *fp;

	maxentry = 1000;
	db = malloc(sizeof(struct database));
	db->dbe = malloc(sizeof(struct dbentry) * maxentry);

	/* Load entries. */
	fp = fopen(dbfile, "r");
	if(fp == NULL){
		perror("load_database");
		return NULL;
	}

	i = 0;
	ch = fgetc(fp);
	do{
		int taxid;

		if(i > maxentry){
			maxentry += 1000;
			db->dbe = realloc(db->dbe, sizeof(struct dbentry) * maxentry);
		}

		ch = read2_integer(fp, ch, &db->dbe[i].taxid);
		if(ch != '\t'){
			fprintf(stderr, "Expect for '\t' after taxon name.\n");
			goto abnormal;
		}

		ch = read_string(fp, &db->dbe[i].name, '\t');
		if(ch != '\t'){
			fprintf(stderr, "Expect for '\t' after taxon name.\n");
			goto abnormal;
		}

		ch = read_string(fp, &db->dbe[i].reffile, '\t');
		db->dbe[i].copy = 1;
		i++;
		ch = fgetc(fp);
	}while(!feof(fp) && ch > 0);

	/* Sort entries by taxonomy id. */
	db->nline = i;
	qsort(db->dbe, db->nline, sizeof(struct dbentry), dbcompar);

	return db;

abnormal:
	if(db){
		if(db->dbe)
			free(db->dbe);
		free(db);
	}

	return NULL;
}

struct dbentry *db_lookup(struct database *db, int taxid)
{
	struct dbentry key, *ret;

	key.taxid = taxid;
	key.name= key.reffile = NULL;
	key.copy = 0;
	ret = bsearch(&key, db->dbe, db->nline, sizeof(struct dbentry), dbcompar);

	return ret;
}

void unload_database(struct database *db)
{
	int i;

	for(i = 0; i < db->nline; i++){
		free(db->dbe[i].name);
		free(db->dbe[i].reffile);
	}

	free(db->dbe);
	free(db);
}

struct abdfile *load_abdfile(char *afile)
{
	struct abdfile *abdf;
	int maxabd, i, ch;
	FILE *fp;

	fp = fopen(afile, "r");
	if(fp == NULL){
		perror("load_adbfile");
		return NULL;
	}

	abdf = malloc(sizeof(struct abdfile));
	memset(abdf, 0, sizeof(struct abdfile));
	maxabd = i = 0;
	ch = fgetc(fp);
	do{
		if(i >= maxabd){
			maxabd += 1000;
			abdf->abd = realloc(abdf->abd, sizeof(struct abund) * maxabd);
		}

		ch = read2_integer(fp, ch, &abdf->abd[i].tax);
		if(ch != '\t'){
			fprintf(stderr, "Expect for '\t' after taxon name.\n");
			goto abnormal;
		}

		ch = read_double(fp, &abdf->abd[i].abund);
		abdf->abd[i].nfrag = 0;
		i++;
		ch = fgetc(fp);
	}while(!feof(fp) && ch > 0);

	abdf->npop = i;
	fclose(fp);

	return abdf;

abnormal:
	if(abdf){
		if(abdf->abd)
			free(abdf->abd);
		free(abdf);
	}
	fclose(fp);

	return NULL;
}

void unload_abdfile(struct abdfile *abdf)
{
	free(abdf->abd);
	free(abdf);
}

int abdcompar(struct abund *a, struct abund *b)
{
	return a->tax - b->tax;
}

int main(int argc, char *argv[])
{
	struct reference *ref;
	struct abdfile *abdf;
	struct abund *abd;
	struct database *db;
	char *afile, *dbfile, *outdir, *prefix, *ptr;
	double *weights, wsum;
	int ntax, npop, nfrags, *nfrag_pop, *nf, i, j, fraglen, rdlen, paired, tax, sd, proflen, nread, *ntrunks;
	char *profpath;

	db = abdf = NULL;
	paired = fraglen = rdlen = 0;
	nread = 1;
	i = 1;
	while(i < argc){
		ptr = argv[i];
		if(*ptr == '-'){
			ptr++;
			switch(*ptr){
				case 'd':	/* Load database. */
					dbfile = argv[++i];
					db = load_database(dbfile);
					break;

				case 'a':	/* Load profile. */
					afile = argv[++i];
					abdf = load_abdfile(afile);
					abd = abdf->abd;
					break;

				case 'o':	/* Output directory. */
					outdir = argv[++i];
					break;

				case 'x':	/* Prefix of output file. */
					prefix = argv[++i];
					break;

				case 'p':	/* Paired-end read. */
					paired = 1;
					nread = 2;
					break;

				case 'f':	/* Read length */
					fraglen = atoi(argv[++i]);
					break;

				case 'r':	/* Read length */
					rdlen = atoi(argv[++i]);
					break;

				case 'n':	/* Total number of reads. */
					nfrags = atoi(argv[++i]);
					break;

				case 't':
					i++;
					break;

				default:
					fprintf(stderr, "Unknown switch -%c\n", *ptr);
					goto abnormal;
			}
		}
		i++;
	}

	if(db == NULL){
		fprintf(stderr, "Reference database is not specified.\n");
		goto abnormal;
	}

	if(abdf == NULL){
		fprintf(stderr, "Abundance profile is not specified.\n");
		goto abnormal;
	}

	if(outdir == NULL){
		fprintf(stderr, "Output directory is not specified.\n");
		goto abnormal;
	}

	if(prefix == NULL){
		fprintf(stderr, "Output prefix is not specified.\n");
		goto abnormal;
	}

	if(fraglen == 0){
		fprintf(stderr, "Fragment length is not specified.\n");
		goto abnormal;
	}

	if((paired && rdlen == 0) || (paired == 0 && rdlen > 0)){
		fprintf(stderr, "-p and -l must be used simultaneously.\n");
		goto abnormal;
	}

	if(paired == 0)
		rdlen = fraglen;
	sd = time(NULL);
//	sd = 1525493755;
//	srand(sd);
	init_rand(sd);
//	fprintf(stderr, "seed=%d\n", sd);

	proflen = strlen(outdir) + strlen(prefix) + 1 + 20 + 4 + 5 + 4 + 5 + 1;	// The constants stands for: '/', taxid, ".tax", chromosome id, ".chr", ".prof", '\0'
	profpath = malloc(sizeof(char) * proflen);

	npop = abdf->npop;
	weights = malloc(sizeof(double) * npop);
	ntrunks = malloc(sizeof(int) * npop);
	memset(ntrunks, 0, sizeof(int) * npop);
	i = 1;
	while(i < argc){
		char *begptr, *endptr;

		ptr = argv[i];
		if(*ptr == '-'){
			ptr++;
			switch(*ptr){
				case 't':
					begptr = argv[++i];
					for(j = 0; j < npop; j++){
						ntrunks[j] = strtol(begptr, &endptr, 10);
						begptr = endptr + 1;
					}

					break;
			}
		}
		i++;
	}

	/* Compute weight of each subpopulation. */
	wsum = 0;
//	fprintf(stderr, "taxid\tlength\tcopy number\tabundance\tweight\n");
	for(i = 0; i < npop; i++){
		struct reference *ref;
		struct dbentry *dbe;
		int len, tax;
		double w;

		tax = abd[i].tax;
		dbe = db_lookup(db, tax);
		ref = load_reference(dbe->reffile);

		len = 0;
		for(j = 0; j < ref->nchr; j++)
			len += ref->chrlen[j];
		w = len * dbe->copy * abd[i].abund;
		weights[i] = (wsum += w);
//		fprintf(stderr, "%d\t%d\t%d\t%.6f\t%.6f\n", tax, len, dbe->copy, abd[i].abund, w);
		unload_reference(ref);
	}

	for(i = 0; i < npop; i++)
		weights[i] /= wsum;

/*	fprintf(stderr, "normalized weights:%.6f", weights[0]);
	for(i = 1; i < npop; i++)
		fprintf(stderr, ", %.6f", weights[i] - weights[i - 1]);
	fprintf(stderr, "\n");*/

	/* Determine source of each read. */
	for(i = 0; i < nfrags; i++){
		double u;

		u = dunif01();
//		u = (double)rand() / RAND_MAX;
		for(j = 0; j < npop; j++)
			if(weights[j] > u)
				break;
		abd[j].nfrag++;
	}

	/* Collect populations of the same taxonomy. */
	qsort(abd, npop, sizeof(struct abund), abdcompar);
	nfrag_pop = malloc(sizeof(int) * npop);
	for(i = 0; i < npop; i++)
		nfrag_pop[i] = abd[i].nfrag;

/*	fprintf(stderr, "New order:%d", abd[0].tax);
	for(i = 1; i < npop; i++)
		fprintf(stderr, ", %d", abd[i].tax);
	fprintf(stderr, "\n");

	fprintf(stderr, "Number of fragments:%d", nfrag_pop[0]);
	for(i = 1; i < npop; i++)
		fprintf(stderr, ", %d", nfrag_pop[i]);
	fprintf(stderr, "\n");*/

	tax = abd[0].tax;
	ntax = 1;
	for(i = 1; i < npop; i++){
		if(tax != abd[i].tax){
			ntax++;
			tax = abd[i].tax;
		}
	}

	/* For each taxon and chromosome, generate a read profile. */
	nf = nfrag_pop;
	for(i = j = 0; i < ntax; i++){
		struct reference *ref;
		struct dbentry *dbe;
		int npop_tax, nchr, *nfrag_cp, *nf2, len, k;
		double *wchr;

		npop_tax = 1;
		tax = abd[j].tax;
		while(++j < npop && abd[j].tax == tax) npop_tax++;

		/* Calculate number of reads of each chromosome */
		dbe = db_lookup(db, tax);
		ref = load_reference(dbe->reffile);
		nchr = ref->nchr;
		wchr = malloc(sizeof(double) * nchr);
		len = 0;
//		fprintf(stderr, "taxon %d: npop=%d, nchr=%d\n", tax, npop_tax, nchr);
		for(k = 0; k < nchr; k++)
			wchr[k] = (len += ref->chrlen[k]);
/*		fprintf(stderr, "chromosome lengths:");
		for(k = 0; k < nchr - 1; k++)
			fprintf(stderr, "%d, ", ref->chrlen[k]);
		fprintf(stderr, "%d\n", ref->chrlen[k]);*/

		for(k = 0; k < nchr; k++)
			wchr[k] /= len;

/*		fprintf(stderr, "chromosome weights:");
		for(k = 0; k < nchr - 1; k++)
			fprintf(stderr, "%.6f, ", wchr[k]);
		fprintf(stderr, "%.6f\n", wchr[k]);*/

		nfrag_cp = malloc(sizeof(int) * nchr * npop_tax);
		memset(nfrag_cp, 0, sizeof(int) * nchr * npop_tax);
		for(k = 0; k < npop_tax; k++){
			int l;

			for(l = 0; l < nf[k]; l++){
				double u;
				int m;

//				u = (double)rand() / RAND_MAX;
				u = dunif01();
				for(m = 0; m < nchr; m++)
					if(wchr[m] > u)
						break;
				nfrag_cp[m * npop_tax + k]++;
			}
		}

/*		for(k = 0; k < nchr; k++){
			int l;

			fprintf(stderr, "fragments of chromosome %d: ", k);
			for(l = 0; l < npop_tax - 1; l++)
				fprintf(stderr, "%d, ", nfrag_cp[k * npop_tax + l]);
			fprintf(stderr, "%d\n", nfrag_cp[k * npop_tax + l]);
		}*/

		/* Load reference of this taxon. */
		nf2 = nfrag_cp;
		for(k = 0; k < nchr; k++){
			struct profile *prof;
			FILE *proffp;

			sprintf(profpath, "%s/%s.tax%d.chr%d.prof", outdir, prefix, tax, k);
			proffp = fopen(profpath, "w+");
			if(proffp == NULL){
				free(wchr);
				free(nfrag_cp);
				unload_reference(ref);
				perror(profpath);
				goto abnormal;
			}

			prof = generate_profile(dbe->reffile, k, npop_tax, nf2, fraglen, paired, rdlen, ntrunks);
			print_profile(prof, proffp);
			unload_profile(prof);
			fclose(proffp);
			nf2 += npop_tax;
		}

		free(wchr);
		free(nfrag_cp);
		unload_reference(ref);
		nf += npop_tax;
	}

	free(nfrag_pop);
	free(ntrunks);
	free(weights);
	unload_abdfile(abdf);
	unload_database(db);
	free(profpath);

	return 0;

abnormal:
	if(nfrag_pop)
		free(nfrag_pop);

	if(weights)
		free(weights);

	if(abdf)
		unload_abdfile(abdf);

	if(db)
		unload_database(db);

	if(profpath)
		free(profpath);

	usage(argv[0]);
	return -1;
}

