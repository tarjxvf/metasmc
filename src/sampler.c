/* The program in this file generate positions and lengths of reads sampled from a reference genome.
 * These information will be passed to coalescent simulator to generate read sequences. */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "global.h"

void usage(char *prog)
{
	fprintf(stderr, "Usage: %s nfrags -g genome_file -f frag_size [-n npop -s size_1 size_2 ... -p -l read_len]\n", prog);
	fprintf(stderr, "nfrags: Number of fragments to be generated.\n");
	fprintf(stderr, "-g genome_file: FASTA file of reference genome.\n");
	fprintf(stderr, "-f frag_size: Fragment size.\n");
	fprintf(stderr, "-n npop: Number of subpopulations. Default value is 1.\n");
	fprintf(stderr, "-s size_1,size_2,...: Relative size of population $2,3...$. Note that the size of population 1 is always 1. If there is npop populations, user should specify npop-1 sizes.\n");
	fprintf(stderr, "-p: If this switch is used, sampler will generate paired-end reads. Default value is 0.\n");
	fprintf(stderr, "-l read_len: If paired-read is switched on, user need to specify read size by this switch. If paired-end is not used, read length will equal to fragment length.\n");
}

int main(int argc, char *argv[])
{
	int nfrag, npop, paired, nread, reflen, fraglen, rdlen, i, j;
	struct frag *fgset;
	char *ptr, *file;
	double *size, allsize;
	unsigned int sd;
	FILE *filp;

	filp = NULL;
	size = NULL;

	if(argc < 2)
		goto abnormal;

	nfrag = atoi(argv[1]);
	if(nfrag <= 0)
		goto abnormal;
	paired = npop = rdlen = 0;
	nread = 1;
	i = 2;
	while(i < argc){
		int ch;

		ptr = argv[i];
		if(*ptr == '-'){
			ptr++;
			switch(*ptr){
				case 's':
					break;

				case 'n':
					npop = atoi(argv[++i]);
					break;

				case 'g':
					file = argv[++i];
					filp = fopen(file, "r");
					if(filp == NULL){
						perror(file);
						goto abnormal;
					}

					/* Read header of first sequence. */
					while((ch = fgetc(filp)) != EOF && ch != '\n');

					/* Count length of reference genome. */
					reflen = 0;
					while((ch = fgetc(filp)) != EOF && ch != '>'){
						if(nucl_index(ch) >= 0)
							reflen++;
					}

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

				default:
					fprintf(stderr, "Unknown switch -%c\n", *ptr);
					goto abnormal;
			}
		}
		i++;
	}

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
	i = 2;
	size = malloc(sizeof(double) * npop);
	size[0] = allsize = 1;
	j = 0;
	while(i < argc){
		char *begptr, *endptr;
		int ch;

		ptr = argv[i];
		if(*ptr == '-'){
			ptr++;
			switch(*ptr){
				case 's':
					begptr = argv[++i];
					for(j = 1; j < npop; j++){
						size[j] = strtod(begptr, &endptr);
						if(size[j] <= 0){
							fprintf(stderr, "Negative population size %6f.\n", size[j]);
							goto abnormal;
						}

						if(j < npop - 1 && *endptr != ','){
							fprintf(stderr, "Expect for ',' between population sizes");
							goto abnormal;
						}

						begptr = endptr + 1;
						allsize += size[j];
					}
					break;
			}
		}
		i++;
	}

	if(npop > 1 && j < npop){
		fprintf(stderr, "%d subpopulation size need to be specified.\n", npop);
		goto abnormal;
	}

	printf("%s\n", file);

	printf("%d,%d\n", reflen, nfrag);
	printf("%d", npop);
	for(i = 1; i < npop; i++)
		printf(",%.2f", size[i]);
	printf("\n");

	sd = time(NULL);
	init_rand(sd);
	for(i = 0; i < nfrag; i++){
		int fragpos, readpos, pop;
		double upop;

		upop = dunif01() * allsize;
		for(pop = 0; pop < npop; pop++){
			if(upop < size[pop])
				break;
			else
				upop -= size[pop];
		}

		fragpos = dunif(reflen - fraglen);
		printf("%d\t%d\t%d\t%d\t%d", i + 1, pop, fragpos, fragpos + fraglen, nread);
		readpos = fragpos;
		for(j = 0; j < nread; j++){
			printf("\t%d\t%d", readpos, readpos + rdlen);
			readpos = fragpos + fraglen - rdlen;
		}
		printf("\n");
	}

	free(size);
	fclose(filp);
	return 0;

abnormal:
	if(size)
		free(size);
	usage(argv[0]);
	return -1;
}
