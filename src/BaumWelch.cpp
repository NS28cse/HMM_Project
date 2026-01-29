#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define	MAXSAMPLE	1000
#define	MAXSTATE	6    // Change 4 -> 6
#define	MAXSYMBOL	10
#define	MAXLEN		512
#define	ERROR		0.01

/*
a[i][j] : transition probability from the state i to j
b[i][k] : output probability for the symbol k at the state i
pi[i] : initial state probability of the state i
c[n][t][i][j] : "gamma" total transition probability from the state i to j
at time t of the sample n
alpha[n][t][i] : forward probability at the state i and time t of the sample n
beta[n][t][i] : backward probability ath teh state i and time t of the sample n
o[n][t] : output symbol at time t of the sample n
*/
int	maxlen[MAXSAMPLE];		/* ?��T?��?��?��v?��?��?��t?��@?��C?��?��?��?��?��Ƃ̒�?��?��?��?��?��i?��[ */
double	a[MAXSTATE][MAXSTATE];
double	na[MAXSTATE][MAXSTATE];		/* a?��̍X?��V?��p */
double	b[MAXSTATE][MAXSYMBOL];
double	nb[MAXSTATE][MAXSYMBOL];	/* b?��̍X?��V?��p */
double	pi[MAXSTATE];
double	npi[MAXSTATE];				/* pi?��̍X?��V?��p */
double	alpha[MAXSAMPLE][MAXLEN][MAXSTATE];
double	beta[MAXSAMPLE][MAXLEN][MAXSTATE];
double	c[MAXSAMPLE][MAXLEN][MAXSTATE][MAXSTATE];
int	o[MAXSAMPLE][MAXLEN];
long totalLen = 0;    // count totalLen for BIC.
int step = 0;    // count step.
int states;

int HMMprint(char *mn)
{
	int	i, j, f;
	FILE	*markovfile;
	char	lognum[MAXLEN];
	char	markovname[MAXLEN];

	strncpy(markovname, mn, MAXLEN);
	if ((markovfile = fopen(markovname, "w")) == NULL)
	{
		printf("markov file write open error : %s\n", markovname);
		exit(0);
	}
	for (i = 0; i < states; i++)
	{
		fprintf(markovfile, "State %c %g\n", i + 'A', pi[i]);
		fprintf(markovfile, "Output");
		for (j = 0; j < MAXSYMBOL; j++)
		{
			fprintf(markovfile, " %g", b[i][j]);
		}
		fprintf(markovfile, "\n");
		fprintf(markovfile, "Transition");
		for (j = 0; j < states; j++)
		{
			fprintf(markovfile, " %g", a[i][j]);
		}
		fprintf(markovfile, "\n");
	}
	fclose(markovfile);

	return 0;
}

void forward_algo()
{
	int		i, j, k, n;

	for (n = 0; n < MAXSAMPLE; n++)
	{
		for (i = 0; i < states; i++)
		{
//			.....
			// Initialize.
			if (maxlen[n] > 0) {
				alpha[n][0][i] = pi[i] * b[i][o[n][0]];
			}
			
		}
	}

	for (n = 0; n < MAXSAMPLE; n++)
	{
		for (k = 0; k + 1 < maxlen[n]; k++)
		{
			for (j = 0; j < states; j++)
			{
//				.....
				// Recursion step for the forward algorithm.
				double sum = 0.0;
				for (i = 0; i < states; i++) {
					sum += alpha[n][k][i] * a[i][j];
				}
				alpha[n][k + 1][j] = sum * b[j][o[n][k + 1]];
			}
		}
	}
}

void backward_algo()
{
	int	i, j, k, n;

	for (n = 0; n < MAXSAMPLE; n++)
	{
		for (i = 0; i < states; i++)
		{
//			.....
			// Initialize.
			if (maxlen[n] > 0) {
				beta[n][maxlen[n] - 1][i] = 1.0;
			}
		}
	}

	for (n = 0; n < MAXSAMPLE; n++)
	{
		for (k = maxlen[n] - 1; k > 0; k--)
		{
			for (i = 0; i < states; i++)
			{
//				.....
				// Induction.
				double sum = 0.0;
				for (j = 0; j < states; j++) {
					sum += a[i][j] * b[j][o[n][k]] * beta[n][k][j];
				}
				beta[n][k - 1][i] = sum;
			}
		}
	}
}

void gamma_algo()
{
	int		i, j, k, l, n;
	double	bunbo;

//	.....
	for (n = 0; n < MAXSAMPLE; n++) {
		if (maxlen[n] <= 1) continue;

		bunbo = 0.0;
		for (i = 0; i < states; i++) {
			bunbo += alpha[n][maxlen[n] - 1][i];
		}

		if (bunbo == 0.0){
            // Clear c if probability is zero to prevent stale data usage.
            for (k = 0; k < maxlen[n] - 1; k++) {
                for (i = 0; i < states; i++) {
                    for (j = 0; j < states; j++) {
                        c[n][k][i][j] = 0.0;
                    }
                }
            }
            continue;
        }

		for (k = 0; k < maxlen[n] - 1; k++) {
			for (i = 0; i < states; i++) {
				for (j = 0; j < states; j++) {
					c[n][k][i][j] = (alpha[n][k][i] * a[i][j] * b[j][o[n][k + 1]] * beta[n][k + 1][j]) / bunbo;
				}
			}
		}
	}

}

void sort_states(int states) {
    int i, j, k, tmp_idx;
    int p[MAXSTATE];
    for (i = 0; i < states; i++) p[i] = i;

    /* Bubble sort based on self-transition probability. */
    for (i = 0; i < states - 1; i++) {
        for (j = i + 1; j < states; j++) {
            if (a[p[i]][p[i]] < a[p[j]][p[j]]) {
                tmp_idx = p[i]; p[i] = p[j]; p[j] = tmp_idx;
            }
        }
    }
    /* After finding new order 'p', rearrange A, B, and PI synchronously. */
    /* Note: Rearrange both rows and columns for matrix A. */
	/* After determining the new order 'p'... */
	/* Swap PI */
	double new_pi[MAXSTATE];
	for(i=0; i<states; i++) new_pi[i] = pi[p[i]];

	/* Swap B */
	double new_b[MAXSTATE][MAXSYMBOL];
	for(i=0; i<states; i++)
    for(k=0; k<MAXSYMBOL; k++) new_b[i][k] = b[p[i]][k];

	/* Swap A (Rows AND Columns) */
	double new_a[MAXSTATE][MAXSTATE];
	for(i=0; i<states; i++)
    for(j=0; j<states; j++) new_a[i][j] = a[p[i]][p[j]];

	for(i=0; i<states; i++) {
    pi[i] = new_pi[i];
    for(k=0; k<MAXSYMBOL; k++) b[i][k] = new_b[i][k];
    for(j=0; j<states; j++) a[i][j] = new_a[i][j];
}
}

double calculate_log_likelihood(int samplecount, int states) {
    double total_loglik = 0.0;
    int n, i;
    for (n = 0; n < samplecount; n++) {
        double prob = 0.0;
        int T = maxlen[n];
        for (i = 0; i < states; i++) {
            prob += alpha[n][T - 1][i]; /* Sum of forward probabilities at time T. */
        }
        total_loglik += log(prob); /* Add log-probability of each sample. */
    }
    return total_loglik;
}

int main(int argc, char *argv[])
{
	char *sampledirname = argv[1];    // Store sample name.
	char *outputdir = argv[2];     // Store output directory path.
	states = atoi(argv[3]);    /* Set number of states from argument. */
	int seed = atoi(argv[4]);      /* Set random seed from argument. */
	char outputfilename[MAXLEN];
	sprintf(outputfilename, "%s/markov_output_s%d_%d.txt", outputdir, states, seed);

	srand((unsigned int)seed);      // Initialize random seed.

	FILE	*samplefile;
//	char	sampledirname[MAXLEN] = "c:/tmp/sample0";
	char	samplefilename[MAXLEN];
//	char	outputfilename[MAXLEN] = "c:/tmp/markov_output.txt";
	char	tmp[MAXLEN];

	int	i, j, k, l, m, n;
	int	readsymbol;
	int	endflag = 0;
	int	converg = 0;
	int	samplecount, timecount;
	double	bunbo;
	double	ibunbo[MAXSTATE];
	double	err;



	for (samplecount = 0, i = 0; i < MAXSAMPLE; i++)
	{


		strncpy(samplefilename, sampledirname, MAXLEN);
		strncat(samplefilename, "/", MAXLEN);
		sprintf(tmp, "%d", i + 1);
		strncat(samplefilename, tmp, MAXLEN);
		strncat(samplefilename, ".txt", MAXLEN);
		// samplefilename == "(sample dir name)\\(i+1).txt"
		if ((samplefile = fopen(samplefilename, "r")) == NULL)
			continue;
		printf("open file : %s\n", samplefilename);
		for (timecount = 0, j = 0; fscanf(samplefile, "%d", &readsymbol) != EOF; j++)
		{
			if (readsymbol >= 0 && readsymbol < MAXSYMBOL)
			{
				o[samplecount][timecount] = readsymbol;
				timecount++;
			}
			else
				printf("symbol range error at sample %d, time %d\n", i, j);
		}
		maxlen[samplecount] = timecount;
        totalLen += timecount; /* Accumulate total length for BIC. */
		samplecount++;
		fclose(samplefile);
	}
	/*
	// Initial setting
	printf("Init a , b and pi\n");
	for (i = 0; i < states; i++)
	{
		for (j = 0; j < states; j++)
		{
			a[i][j] = 1 / (double)MAXSTATE;
		}
		pi[i] = 1 / (double)MAXSTATE;
		for (j = 0, bunbo = 0; j < MAXSYMBOL; j++)
		{
			b[i][j] = (double)((j + 1) % (i + 1) + 1);
			bunbo += b[i][j];
		}
		for (j = 0; j < MAXSYMBOL; j++)
			b[i][j] /= bunbo;
	}
	*/
	// Initial setting.
	printf("Init a , b and pi with random values\n");

	// Initialize pi with random values and normalize.
	bunbo = 0;
	for (i = 0; i < states; i++) {
		pi[i] = (double)rand() / RAND_MAX;
		bunbo += pi[i];
	}
	for (i = 0; i < states; i++) pi[i] /= bunbo;

	// Initialize a (Transition probability) with random values and normalize.
	for (i = 0; i < states; i++)
	{
		bunbo = 0;
		for (j = 0; j < states; j++)
		{
			a[i][j] = (double)rand() / RAND_MAX;
			bunbo += a[i][j];
		}
		for (j = 0; j < states; j++) a[i][j] /= bunbo;

		// Initialize b (Output probability) with random values and normalize.
		bunbo = 0;
		for (j = 0; j < MAXSYMBOL; j++)
		{
			b[i][j] = (double)rand() / RAND_MAX;
			bunbo += b[i][j];
		}
		for (j = 0; j < MAXSYMBOL; j++) b[i][j] /= bunbo;
	}
	/* Perform the first forward-backward pass to get initial alpha/beta values for Step 0. */
	forward_algo();
	double old_loglik = calculate_log_likelihood(samplecount, states); 
	printf("LABEL Step LogLik DError\n"); /* Print the header for history output. */

	// loop start
	converg = 0;
	while (converg != 1)
	{
		//		printf("forward_algo\n");
		forward_algo();

		//		printf("backward_algo\n");
		backward_algo();

		//		printf("gamma_algo\n");
		gamma_algo();

		//		printf("parameter update : ");
		for (i = 0; i < states; i++)
		{
			// update a
			for (j = 0; j < states; j++)
			{
				na[i][j] = 0;
				bunbo = 0;
				for (n = 0; n < MAXSAMPLE; n++)
				{
					for (k = 0; k < maxlen[n] - 1; k++)
					{
						na[i][j] += c[n][k][i][j];
						for (l = 0; l < states; l++)
							bunbo += c[n][k][i][l];
					}
				}
				if (bunbo == 0)
				{
					if (na[i][j] != 0)
					{
						printf("bunbo a=0:%d,%d : ", i, j);
						na[i][j] = 0;
					}
				}
				else
				{
					na[i][j] /= bunbo;
				}
			}

			// update pi
			npi[i] = 0;
			for (n = 0; n < MAXSAMPLE; n++)
			{
				for (j = 0; j < states; j++)
					npi[i] += c[n][0][i][j];
			}

			// update b
			for (j = 0; j < MAXSYMBOL; j++)
			{
				nb[i][j] = 0;
				bunbo = 0;
				for (n = 0; n < MAXSAMPLE; n++)
				{
					/*
					for (k = 0; k < maxlen[n]; k++)
					{
						if (o[n][k] == j)
							for (l = 0; l < states; l++)
								nb[i][j] += c[n][k][i][l];
						for (l = 0; l < states; l++)
							bunbo += c[n][k][i][l];
					}
					*/
					for (k = 0; n < maxlen[n]; k++)
					{
						double gamma_val = 0.0;
						// Calculation of gamma at time k
						if (k < maxlen[n] - 1)
						{
							// For t < T, gamma = sum(xi)
							for (l = 0; l < states; l++)
								gamma_val += c[n][k][i][l];
						}
						else
						{
							// For t = T, gamma = alpha / P(O)
							// Re-calculate P(O) (bunbo in gamma_algo) for this sample
							double p_O = 0.0;
							for (l = 0; l < states; l++)
								p_O += alpha[n][k][l];
							if (p_O != 0)
								gamma_val = alpha[n][k][i] / p_O;
						}
						if (o[n][k] == j)
						    nb[i][j] += gamma_val;
						bunbo += gamma_val;
					}
				}
				if (bunbo == 0)
				{
					if (nb[i][j] != 0)
						nb[i][j] = 0;
				}
				else
					nb[i][j] /= bunbo;
			}
		}
		// update pi
		for (bunbo = 0, i = 0; i < states; i++)
		{
			bunbo += npi[i];
		}
		for (i = 0; i < states; i++)
			npi[i] /= bunbo;
		//	printf("continue -> ");
		for (i = 0; i < states; i++)
		{
			for (j = 0; j < states; j++)
			{
				err = fabs(a[i][j] - na[i][j]);
				if (err / a[i][j] > ERROR)
				{
					printf("error on a[%d][%d] = %g %g > %g\n", i, j, err, err / a[i][j], ERROR);
					break;
				}
			}
			if (j < states)
			{
				break;
			}

			for (j = 0; j < MAXSYMBOL; j++)
			{
				err = fabs(b[i][j] - nb[i][j]);
				if (err / b[i][j] > ERROR)
				{
					printf("error on b[%d][%d] = %g %g > %g\n", i, j, err, err / b[i][j], ERROR);
					break;
				}
			}
			if (j < MAXSYMBOL)
			{
				break;
			}

			err = fabs(pi[i] - npi[i]);
			if (err / pi[i] > ERROR)
			{
				printf("error on pi[%d] = %g %g > %g\n", i, err, err / pi[i], ERROR);
				break;
			}
		}
		if (!(i < states))
		{
			converg = 1;
		}
		else
		{
			converg = 0;
		}
		for (i = 0; i < states; i++)
		{
			for (j = 0; j < states; j++)
			{
				a[i][j] = na[i][j];
			}
			for (j = 0; j < MAXSYMBOL; j++)
			{
				b[i][j] = nb[i][j];
			}
			pi[i] = npi[i];
		}
		step++; /* Increment step count. */
		forward_algo(); /* Prepare alpha for likelihood calculation. */
		double current_loglik = calculate_log_likelihood(samplecount, states);
		printf("HISTORY %d %f %f\n", step, current_loglik, current_loglik - old_loglik); /* Print history. */
		old_loglik = current_loglik; /* Update likelihood for next step. */
	}
	printf("loop end and print HMM\n");
	sort_states(states); /* Finalize state order. */

	printf("LABEL SampleName States Seed Symbols TotalLen FinalLogLik Iterations\n");
	printf("RESULTS %s %d %d %d %ld %f %d\n", 
        sampledirname, states, seed, MAXSYMBOL, totalLen, old_loglik, step); /* Output final result for MATLAB. */
    
    char finalPath[MAXLEN];
    sprintf(finalPath, "%s/markov_output_s%d_%d.txt", outputdir, states, seed); /* Generate file path. */
	HMMprint(outputfilename);
}
