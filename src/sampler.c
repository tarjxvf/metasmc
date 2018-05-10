/* The program in this file generate positions and lengths of reads sampled from a reference genome.
 * These information will be passed to coalescent simulator to generate read sequences. */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include "global.h"

void usage(char *prog)
{
	fprintf(stderr, "Usage: %s -g genome_file -c chr_num -f frag_size -r nfrag_1,nfrag_2,... [-n npop -p -l read_len -t ntrunk1,ntrunk2,...]\n", prog);
	fprintf(stderr, "-g genome_file: FASTA file of reference genome.\n");
	fprintf(stderr, "-c chr_num: 0-based index of chromosome in genome_file. Default is 0.\n");
	fprintf(stderr, "-f frag_size: Fragment size.\n");
	fprintf(stderr, "-n npop: Number of subpopulations. Default value is 1.\n");
	fprintf(stderr, "-r nfrag_1,nfrag_2,...: Number of reads of each subpopulations.\n");
	fprintf(stderr, "-p: If this switch is used, sampler will generate paired-end reads. Default value is 0.\n");
	fprintf(stderr, "-l read_len: If paired-read is switched on, user need to specify read size by this switch. If paired-end is not used, read length will equal to fragment length.\n");
	fprintf(stderr, "-t ntrunk1,ntrunk2,...: Number of whole chromosome in each subpopulation.");
}

int main(int argc, char *argv[])
{
	int nfrag, npop, paired, nread, reflen, fraglen, rdlen, i, j, *nfrag_pop, chrnum, *ntrunks;
	struct reference *ref;
	struct profile *prof;
	struct frag *fgset;
	char *ptr, *file;
	double *size, allsize;
	unsigned int sd;
	FILE *filp;

	filp = NULL;
	nfrag_pop = NULL;
	ref = NULL;
//	size = NULL;

	if(argc < 2)
		goto abnormal;

//	nfrag = atoi(argv[1]);
//	if(nfrag <= 0)
//		goto abnormal;
	paired = npop = rdlen = 0;
	nread = 1;
	i = 1;
	chrnum = 0;
	while(i < argc){
		int ch;

		ptr = argv[i];
		if(*ptr == '-'){
			ptr++;
			switch(*ptr){
				case 'r':
					i++;
					break;

				case 'c':
					chrnum = atoi(argv[++i]);
					break;

				case 'n':
					npop = atoi(argv[++i]);
					break;

				case 'g':
					file = argv[++i];
//					filp = fopen(file, "r");
//					if(filp == NULL){
//						perror(file);
//						goto abnormal;
//					}

					ref = load_reference(file);

					/* Read header of first sequence. */
//					while((ch = fgetc(filp)) != EOF && ch != '\n');

					/* Count length of reference genome. */
//					reflen = 0;
//					while((ch = fgetc(filp)) != EOF && ch != '>'){
//						if(nucl_index(ch) >= 0)
//							reflen++;
//					}

					break;

				case 'f':
					fraglen = atoi(argv[++i]);
					break;

				case 'p':
					paired = 1;
					nread = 2;
					break;

				case 'l':
					rdlen = atoi(argv[++i]);
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

	if(ref == NULL){
		fprintf(stderr, "Reference file is not specified.\n");
		goto abnormal;
	}

	if(chrnum > ref->nchr){
		fprintf(stderr, "Chromosome index is too large.\n");
		goto abnormal;
	}
	reflen = ref->chrlen[chrnum];

	if(file == NULL)
		goto abnormal;

	if((paired && rdlen == 0) || (paired == 0 && rdlen > 0)){
		fprintf(stderr, "-p and -l must be used simultaneously.\n");
		goto abnormal;
	}

	if(fraglen == 0)
		goto abnormal;

	if(rdlen == 0)
		rdlen = fraglen;

	if(npop == 0)
		npop = 1;
	i = 1;
//	size = malloc(sizeof(double) * npop);
//	size[0] = allsize = 1;
	nfrag_pop = malloc(sizeof(int) * npop);
	ntrunks = malloc(sizeof(int) * npop);
	memset(ntrunks, 0, sizeof(int) * npop);
	j = 0;
	while(i < argc){
		char *begptr, *endptr;
		int ch;

		ptr = argv[i];
		if(*ptr == '-'){
			ptr++;
			switch(*ptr){
				case 'r':
					begptr = argv[++i];
					nfrag = 0;
					for(j = 0; j < npop; j++){
						nfrag_pop[j] = strtol(begptr, &endptr, 10);
						nfrag += nfrag_pop[j];
						begptr = endptr + 1;
					}
					break;

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
	unload_reference(ref);

	if(npop > 1 && j < npop){
		fprintf(stderr, "%d subpopulation size need to be specified.\n", npop);
		goto abnormal;
	}

	sd = time(NULL);
	init_rand(sd);
	prof = generate_profile(file, chrnum, npop, nfrag_pop, fraglen, paired, rdlen, ntrunks);
	print_profile(prof, stdout);
	unload_profile(prof);

//	free(size);
	free(nfrag_pop);
//	fclose(filp);
	return 0;

abnormal:
	if(nfrag_pop)
		free(nfrag_pop);
//	if(size)
//		free(size);
	usage(argv[0]);
	return -1;
}
