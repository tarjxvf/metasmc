/* This is the main program of coalescent simulator. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

#include "global.h"
#include "smc.h"
#include "mutation.h"
#include "slab.h"

#define NEXT_NOBLANK(fp, ch)    while(isblank((ch) = fgetc(fp)))

/*int read_integer(FILE *filp, int *val)
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

	return ch;
}*/

void usage(char *prog)
{
	fprintf(stderr, "Usage: %s -i read_profile -t theta [Other options]\n", prog);
	fprintf(stderr, "-i read_profile : read profile consisting of path to reference genome, number of fragments and position of each fragment. If this switch is not used, the input are read from stdandard input\n");
	fprintf(stderr, "-s seed : Seed of pseudorandom number generator.\n");
	fprintf(stderr, "-F maxfrag : Maximum number of clones per reigion.\n");
	fprintf(stderr, "-t theta1,theta2,... : per-locus mutation rate\n");
	fprintf(stderr, "-r rho : per-locus recombination rate\n");
	fprintf(stderr, "-T tree_file: Print coalescent trees in Newick format to tree_file.\n");
	fprintf(stderr, "-o read_file: Print reads in FASTA format to read_file.\n");
	fprintf(stderr, "-d tdiv : divergence time from reference genome.\n");
	fprintf(stderr, "-m model1,model2,...: use DNA evolution models. Available models are JC69, K80, F81, HKY85, T92, TN93 and GTR. The description of mutation models can be found in http://en.wikipedia.org/wiki/Models_of_DNA_evolution. The default mutation model is infinite-site model.\n");
	fprintf(stderr, "-f freq1_A,freq1_C,freq1_G,freq1_T;freq2_A,freq2_C,freq2_G,freq2_T,...: Equalibrium allele frequencies.\n");
	fprintf(stderr, "-R pars1;pars2;...: Set parameters of mutation model. Number of parameters depends on mutation model.\n");
	fprintf(stderr, "-G alpha: Default exponential growth rate (at time 0)\n");
	fprintf(stderr, "-g i alpha: Default exponential growth rate (at time 0) of subpopulation i\n");
	fprintf(stderr, "-S size_1,size_2,...: Size of each subpopulation at time 0. All sizes are set to 1 if this switch is ignored.\n");
	fprintf(stderr, "-M m_12,m_13,...,m_21,m_23,...: Default migration matrix (at time 0)\n");
	fprintf(stderr, "-eG t alpha  (Modify growth rate of all pop's.)\n");
	fprintf(stderr, "-eg t i alpha_i  (Modify growth rate of pop i.)\n");
	fprintf(stderr, "-eM t mig_rate   (Modify the mig matrix so all elements are mig_rate/(npop-1)\n");
	fprintf(stderr, "-em t i j m_ij    (i,j-th element of mig matrix set to m_ij at time t)\n");
	fprintf(stderr, "-eN t size  (Modify pop sizes. New sizes = size*N0)\n");
	fprintf(stderr, "-en t i size_i  (Modify pop size of pop i.  New size of popi = size_i*N0 .)\n");
	fprintf(stderr, "-es t i proportion  (Split: pop i -> pop-i + pop-npop, npop increases by 1.\n");
	fprintf(stderr, "        proportion is probability that each lineage stays in pop-i. (p, 1-p are admixt. proport.\n");
	fprintf(stderr, "        Size of pop npop is set to N0 and alpha = 0.0 , size and alpha of pop i are unchanged.\n");
	fprintf(stderr, "-ej t i j   ( Join lineages in pop i and pop j into pop j\n");
	fprintf(stderr, "        size, alpha and M are unchanged.\n");
}

int main(int argc, char *argv[])
{
	int i, ch, nfrag, pop, npop, nsplt;	/* npar: number of parameters in mutation model,
						 * nsplt: number of subpopulations created by splitting events.
						 * nfrag: number of fragments.
						 * npop: number of subpopulations at time 0. */
	int print_tree;	// Equals to 1 if user want to print local genealogies.
	int maxfrag; // Maximum number of reads in a group
	struct list_head *evl, *l;
	struct genealogy *G;
	struct frag *fgset;
	struct profile *prof;
	struct mutation *mmut;
	struct config *cfg;
	struct event *ev;
	unsigned int sd;
	double rho, tdiv;
	double *grate, *mmig;
	char file[1000], *treefile;
	FILE *filp, *treefp, *readfp;

	/* Parse parameters. */
	filp = treefp = readfp = NULL;	// This will be replaced by stdin if -T switch is used.
	generate_sequence = NULL;
	G = NULL;
	mmut = NULL;
	cfg = 0;
	i = 1;
	nsplt = 0;
	maxfrag = 0;
	grate = NULL;
	mmig = NULL;
	rho = 0;
	prof = NULL;

	print_tree = 0;		// This will be overwritten if -T is used.
	sd = time(NULL);	// By default, seed is generated using time(NULL). This will be overwritten if -s switch is used.
	generate_sequence = generate_sequence_infinite_fast;	// This will be overwritten if -m switch is used.

	// First round: set up basic parameters
	while(i < argc){
		char *ptr, *optarg;
		int j;

		ptr = argv[i];
		if(*ptr == '-'){
			ptr++;
			switch(*ptr){
				case 'i':
					filp = fopen(argv[++i], "r");
					if(filp == NULL){
						perror(argv[i]);
						goto abnormal;
					}

					break;

				case 'r':
					rho = atof(argv[++i]);
					break;

				case 'd':
					tdiv = atof(argv[++i]);
					break;

				case 't': case 'm':	/* Ignore this parameter for the moment. This will be parsed in second round. */
					break;

				case 'R':case 'f':
					/* Ignore this parameter in first run because mutation model may be unknown. */
					break;

				case 'F':
					maxfrag = atoi(argv[++i]);
					break;

				case 'S':
					i++;
					break;

				case 'M':case 'g':case 'G':
				case 'e':
					/* Ignore this parameter in first run because number of subpopulations may be unknown. */
					break;

				case 'T':
					print_tree = 1;
					i++;
					if(i >= argc){
						fprintf(stderr, "Treefile is unspecified.\n");
						goto abnormal;
					}

					treefile = argv[i];
					if(treefile[0] == '-' && treefile[1] == '\0'){
						treefp = stdout;

					}else{
						treefp = fopen(treefile, "w+");
					}

					if(treefp == NULL){
						perror(treefile);
						goto abnormal;
					}

					break;

				case 's':
					sd = atoi(argv[++i]);
					break;

				default:
					if(!isdigit(*ptr)){
						fprintf(stderr, "Unknown switch -%c\n", *ptr);
						goto abnormal;
					}
			}
		}
		i++;
	}

	/* If input file is not assigned, read profile from standard input. */
	if(filp == NULL)
		filp = stdin;

	prof = load_profile(filp);
	npop = prof->npop;
	nfrag = prof->nfrag;

	if(maxfrag == 0)
		maxfrag = nfrag;

	cfg = create_config(sd, print_tree, 1, treefp, readfp, maxfrag, prof, rho);

	// Third round: set up splitting events
	i = 1;
	while(i < argc){
		char *ptr, *optarg, *endptr;
		double t;	// Event time
		int j, k;

		ptr = argv[i];
		if(*ptr == '-'){
			ptr++;

			switch(*ptr){
				case 'e':	/* User-specified demographic event. */
					t = atof(argv[++i]);
					ptr++;
					switch(*ptr){
						double prop;
						int pop;

						case 's':	/* Split a population. */
							pop = strtol(argv[++i], &endptr, 10);
							if(endptr == argv[i]){
								fprintf(stderr, "Invalid parameter for -es\n");
								goto abnormal;
							}

							prop = strtod(argv[++i], &endptr);
							if(endptr == argv[i]){
								fprintf(stderr, "Invalid parameter for -es\n");
								goto abnormal;
							}

							add_event_splt(cfg, t, pop, prop);

							break;
					}

					break;
			}
		}
		i++;
	}

	nsplt = cfg->nsplt;

	mmut = malloc(sizeof(struct mutation) * cfg->npop_all);
	memset(mmut, 0, sizeof(struct mutation) * cfg->npop_all);

	// Second round: set up basic parameters of mutation models
	i = 1;
	while(i < argc){
		char *ptr, *ptr2, *optarg, *begptr, *endptr;
		int j;

		ptr = argv[i];
		if(*ptr == '-'){
			ptr++;
			switch(*ptr){
				case 'i': case 'r': case 'R': case 'f': case 'F': case 'M': case 'G': case 'g': case 'e': case 'T': case 's': case 'S':
					i++;
					break;

				case 't':
					begptr = argv[++i];
					for(j = 0; j < cfg->npop_all; j++){
						mmut[j].theta = strtod(begptr, &endptr);

						if(j < npop - 1 && *endptr != ','){
							break;
//							fprintf(stderr, "Expect for ',' between thetas\n");
//							goto abnormal;
						}
						begptr = endptr + 1;
					}

					if(j > 0 && j < npop){
						fprintf(stderr, "Number of theta parameters should be equal to either 1 or number of populations.\n");
						goto abnormal;

					}else if(j == 0){	// All populations has the same mutation rates at time 0
						for(j = 1; j < cfg->npop_all; j++){
							mmut[j].theta = mmut[0].theta;
						}
					}

					break;

				case 'm':
					/* Detect mutation model. Note that mutation is not assigned here because mutation rate may be unknown */
					ptr2 = argv[++i];
					for(j = 0; j < cfg->npop_all; j++){
						int inc;

						if(!strncmp("JC69", ptr2, 4)){
							mmut[j].model = MODEL_JC69;
							inc = 4;

						}else if(!strncmp("K80", ptr2, 3)){
							mmut[j].model = MODEL_K80;
							inc = 3;

						}else if(!strncmp("F81", ptr2, 3)){
							mmut[j].model = MODEL_F81;
							inc = 3;

						}else if(!strncmp("HKY85", ptr2, 5)){
							mmut[j].model = MODEL_HKY85;
							inc = 5;

						}else if(!strncmp("T92", ptr2, 3)){
							mmut[j].model = MODEL_T92;
							inc = 3;

						}else if(!strncmp("TN93", ptr2, 4)){
							mmut[j].model = MODEL_TN93;
							inc = 4;

						}else if(!strncmp("GTR", ptr2, 3)){
							mmut[j].model = MODEL_GTR;
							inc = 3;

						}else{
							ptr2[inc] = '\0';
							fprintf(stderr, "Unknown mutation model %s\n", ptr2);
							goto abnormal;
						}

						if(j < npop - 1 && ptr2[inc] != ','){
							break;
//							fprintf(stderr, "Expect for ',' between mutation models");
//							goto abnormal;
						}

						ptr2 += inc + 1;
					}

					if(j > 0 && j < npop){
						fprintf(stderr, "Number of models after -m should be equal to either 1 or number of populations.\n");
						goto abnormal;

					}else if(j == 0){
						for(j = 1;j < npop + nsplt; j++){
							mmut[j].model = mmut[0].model;
						}
					}

					generate_sequence = generate_sequence_model_fast;
					break;

				default:
					if(!isdigit(*ptr)){
						fprintf(stderr, "Unknown switch -%c\n", *ptr);
						goto abnormal;
					}
			}
		}
		i++;
	}

	// Set default equilibrium nuclieotide frequencies
	for(i = 0; i < cfg->npop_all; i++){
		int j;

		for(j = 0; j < NUM_NUCS; j++)
			mmut[i].pi[j] = (double)1 / NUM_NUCS;
	}

	fgset = prof->fgset;
	qsort((void *)fgset, (size_t)nfrag, sizeof(struct frag), fgcompar);

	// Fourth round: set up demographic model and additional parameters of mutation models
	grate = malloc(sizeof(double) * npop);
	memset(grate, 0, sizeof(double) * npop);
	mmig = malloc(sizeof(double) * npop * npop);
	memset(mmig, 0, sizeof(double) * npop * npop);
	for(i = 0; i < npop * npop; i++){
		if(i % (npop + 1))
			mmig[i] = (double)1 / (npop - 1);
		else
			mmig[i] = 0;
	}

	i = 1;
	while(i < argc){
		char *ptr, *optarg, *begptr, *endptr;
		double t, *size;	// Event time
		double ftotal;
		int j, k;

		ptr = argv[i];
		if(*ptr == '-'){
			ptr++;

			switch(*ptr){
				case 'i': case 't': case 'r': case 'm': case 'T': case 'F': case 's':
					i++;
					break;

				case 'R':
					begptr = endptr = argv[++i];
					for(j = 0; j < cfg->npop_all; j++){
						for(k = 0; k < npar[mmut[j].model]; k++){
							mmut[j].mpar[k] = strtod(begptr, &endptr);

							if(k < npar[mmut[j].model] - 1 && *endptr != ','){
								fprintf(stderr, "Expect for ',' between model parameters.\n");
								goto abnormal;
							}

							begptr = endptr + 1;
						}

						if(j < npop - 1 && *endptr != ':'){
							break;

						}else{
							begptr = endptr + 1;
						}
					}

					if(j > 0 && j < npop + nsplt){
						fprintf(stderr, "Number of parameter list should be equal to either 1 or number of populations.\n");
						goto abnormal;

					}else if(j == 0){
						for(j = 1; j < cfg->npop_all; j++){
							if(mmut[j].model != mmut[j].model){
								goto abnormal;

							}else{
								memcpy(mmut[j].mpar, mmut[0].mpar, sizeof(double) * npar[mmut[j].model]);
							}
						}
					}

					break;

				case 'f':	/* Load base frequencies. */
					begptr = endptr = argv[++i];
					for(j = 0; j < cfg->npop_all; j++){
						ftotal = 0;
						for(k = 0; k < NUM_NUCS; k++){
							mmut[j].pi[k] = strtod(begptr, &endptr);
							ftotal += mmut[j].pi[k];

							if(k < NUM_NUCS - 1 && *endptr != ','){
								fprintf(stderr, "Expect for ',' between allele frequencies of a population.\n");
								goto abnormal;
							}
							begptr = endptr + 1;
						}

						if(ftotal != 1){
							fprintf(stderr, "Allele frequencies must sum up to 1.\n");
							goto abnormal;
						}

						if(j < npop - 1 && *endptr != ':'){
							break;

						}else{
							begptr = endptr + 1;
						}
					}

					if(j > 0 && j < npop){
						fprintf(stderr, "Number of allele frequency list should be equal to either 1 or number of populations.\n");
						goto abnormal;

					}else if(j == 0){
						for(j = 1; j < cfg->npop_all; j++){
							memcpy(mmut[j].pi, mmut[0].pi, sizeof(double) * NUM_NUCS);
						}
					}

					break;

				case 'S':	/* Specify subpopulation sizes at time 0. */
					begptr = endptr = argv[++i];
					for(j = 0; j < npop; j++){
						cfg->size[j] = strtod(begptr, &endptr);
						begptr = endptr + 1;
					}

					break;

				case 'M':	/* Set default migration matrix. */
					if(npop <= 1){
						fprintf(stderr, "-M must be used in the case of structured population.\n");
						goto abnormal;
					}

					begptr = argv[++i];
					/* Load non-diagonal terms of migration matrix at time 0 */
					for(j = 0; j < npop; j++){
						for(k = 0; k < npop; k++){
							if(j == k){
								mmig[k * prof->npop + k] = 0;	// Diagonal terms are always 0

							}else{
								mmig[j * prof->npop + k] = strtod(begptr, &endptr);
								if(j * npop + k < npop * (npop - 1) && *endptr != ','){
									fprintf(stderr, "Expect for ',' between migration rates.\n");
									goto abnormal;
								}
								begptr = endptr + 1;
							}
						}
					}

					break;

				case 'G':
					grate[0] = strtod(argv[++i], &endptr);
					for(j = 1; j < cfg->npop; j++)
						grate[j] = grate[0];

					set_growth_rates(cfg, grate);

					break;

				case 'g':
					if(npop <= 1){
						fprintf(stderr, "-g can only be used with subdivision.\n");
						goto abnormal;
					}

					j = strtol(argv[++i], &endptr, 10);
					grate[j] = strtod(argv[++i], &endptr);

					break;

				case 'e':	/* User-specified demographic event. */
					t = atof(argv[++i]);
					ptr++;
					switch(*ptr){
						int pop, popi, popj;
						double alpha, prop, size, rmig;

						case 'G':	/* Modify growth rate of all subpupulations. */
							alpha = strtod(argv[++i], &endptr);
							if(endptr == argv[i]){
								fprintf(stderr, "Invalid parameter for -eG\n");
								goto abnormal;
							}
							add_event_ggro(cfg, t, alpha);

							break;

						case 'g':	/* Modify growth rate of a subpopulation. */
							pop = strtol(argv[++i], &endptr, 10);
							if(endptr == argv[i]){
								fprintf(stderr, "Invalid parameter for -eg\n");
								goto abnormal;
							}

							alpha = strtod(argv[++i], &endptr);
							if(endptr == argv[i]){
								fprintf(stderr, "Invalid parameter for -eg\n");
								goto abnormal;
							}
							add_event_grow(cfg, t, pop, alpha);

							break;

						case 'M':	/* Modify migration rate of all populations. */
							rmig = strtod(argv[++i], &endptr);
							if(endptr == argv[i]){
								fprintf(stderr, "Invalid parameter for -eM\n");
								goto abnormal;
							}
							add_event_gmig(cfg, t, rmig);

							break;

						case 'm':	/* Modify migration rate from subpopulation i to subpopulation j. */
							popi = strtol(argv[++i], &endptr, 10);
							if(endptr == argv[i]){
								fprintf(stderr, "Invalid parameter for -em\n");
								goto abnormal;
							}

							popj = strtol(argv[++i], &endptr, 10);
							if(endptr == argv[i]){
								fprintf(stderr, "Invalid parameter for -em\n");
								goto abnormal;
							}

							rmig = strtod(argv[++i], &endptr);
							if(endptr == argv[i]){
								fprintf(stderr, "Invalid parameter for -em\n");
								goto abnormal;
							}
							add_event_rmig(cfg, t, popi, popj, rmig);

							break;

						case 'N':	/* Modify size of all population. */
							size = strtod(argv[++i], &endptr);
							if(endptr == argv[i]){
								fprintf(stderr, "Invalid parameter for -eN\n");
								goto abnormal;
							}
							add_event_gsiz(cfg, t, size);

							break;

						case 'j':	/* Merge subpopulation i into subpopulation j. */
							popi = strtol(argv[++i], &endptr, 10);
							if(endptr == argv[i]){
								fprintf(stderr, "Invalid parameter for -ej\n");
								goto abnormal;
							}

							popj = strtol(argv[++i], &endptr, 10);
							if(endptr == argv[i]){
								fprintf(stderr, "Invalid parameter for -ej\n");
								goto abnormal;
							}
							add_event_join(cfg, t, popi, popj);

							break;

						case 's':	/* Split a population. */
							i += 2;
							break;

						case 'n':	/* Modify size of a population. */
							pop = strtol(argv[++i], &endptr, 10);
							if(endptr == argv[i]){
								fprintf(stderr, "Invalid parameter for -en\n");
								goto abnormal;
							}

							size = strtod(argv[++i], &endptr);
							if(endptr == argv[i]){
								fprintf(stderr, "Invalid parameter for -en\n");
								goto abnormal;
							}
							add_event_size(cfg, t, pop, size);

							break;
					}

					break;

				default:
					fprintf(stderr, "Unknown switch -%c\n", *ptr);
					goto abnormal;
			}
		}
		i++;
	}

	set_growth_rates(cfg, grate);
	free(grate);
	set_migration_matrix(cfg, mmig);
	free(mmig);
	grate = mmig = NULL;

	for(i = 0; i < npop + nsplt; i++){
		register_mutation_model(cfg, i, &mmut[i]);
		init_mutation_model(&mmut[i]);
	}

	G = alloc_genealogy(cfg, prof);

	/* Sort fragments according to start positions. */
	qsort((void *)fgset, (size_t)nfrag, sizeof(struct frag), fgcompar);

	fprintf(stderr, "seed=%d\n", sd);
	init_rand(sd);

//	dump_config(cfg);

	simulate(G, prof);

	/* Print reads */
	for(i = 0; i < prof->nfrag; i++){
		print_fragment(stdout, &prof->fgset[i]);
	}

	free(mmut);

	/* Finalize genealogy */
	clear_genealogy(G);
	destroy_genealogy(G);
	destroy_config(cfg);
	unload_profile(prof);

	fclose(filp);

	return 0;

abnormal:
	if(grate)
		free(grate);

	if(mmig)
		free(mmig);

	if(mmut)
		free(mmut);

	if(treefp)
		fclose(treefp);

	if(filp)
		fclose(filp);

	if(G){
		destroy_genealogy(G);
		free(G);
	}

	if(cfg)
		destroy_config(cfg);

	if(prof)
		unload_profile(prof);

	usage(argv[0]);
	return -1;
}

