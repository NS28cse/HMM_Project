#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>    // for srand().

#define	MAXSAMPLE	1000
#define	MAXSTATE	6     // 4 -> 6.
#define	MAXSYMBOL	10
#define	MAXLEN		512
#define	ERROR		0.01
#define MIN_PROB	1e-100    // small value to avoid divison by zero.

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
int	maxlen[MAXSAMPLE];		/* サンプルファイルごとの長さを格納 */
double	a[MAXSTATE][MAXSTATE];
double	na[MAXSTATE][MAXSTATE];		/* aの更新用 */
double	b[MAXSTATE][MAXSYMBOL];
double	nb[MAXSTATE][MAXSYMBOL];	/* bの更新用 */
double	pi[MAXSTATE];
double	npi[MAXSTATE];				/* piの更新用 */
double	alpha[MAXSAMPLE][MAXLEN][MAXSTATE];
double	beta[MAXSAMPLE][MAXLEN][MAXSTATE];
double	c[MAXSAMPLE][MAXLEN][MAXSTATE][MAXSTATE];
int	o[MAXSAMPLE][MAXLEN];
double scale[MAXSAMPLE][MAXLEN];    // scale value.
double total_log_lik;    
double gamma_val[MAXSAMPLE][MAXLEN][MAXSTATE];

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
	for (i = 0; i < MAXSTATE; i++)
	{
		fprintf(markovfile, "State %c %g\n", i + 'A', pi[i]);
		fprintf(markovfile, "Output");
		for (j = 0; j < MAXSYMBOL; j++)
		{
			fprintf(markovfile, " %g", b[i][j]);
		}
		fprintf(markovfile, "\n");
		fprintf(markovfile, "Transition");
		for (j = 0; j < MAXSTATE; j++)
		{
			fprintf(markovfile, " %g", a[i][j]);
		}
		fprintf(markovfile, "\n");
	}
	fclose(markovfile);
}

void forward_algo()
{
	int		i, j, k, n;

	for (n = 0; n < MAXSAMPLE; n++)
	{
		for (i = 0; i < MAXSTATE; i++)
		{
//			.....
		}
	}

	for (n = 0; n < MAXSAMPLE; n++)
	{
		for (k = 0; k + 1 < maxlen[n]; k++)
		{
			for (j = 0; j < MAXSTATE; j++)
			{
//				.....
			}
		}
	}
}

void backward_algo()
{
	int	i, j, k, n;

	for (n = 0; n < MAXSAMPLE; n++)
	{
		for (i = 0; i < MAXSTATE; i++)
		{
//			.....
		}
	}

	for (n = 0; n < MAXSAMPLE; n++)
	{
		for (k = maxlen[n] - 1; k > 0; k--)
		{
			for (i = 0; i < MAXSTATE; i++)
			{
//				.....
			}
		}
	}
}

void gamma_algo()
{
	int		i, j, k, l, n;
	double	bunbo;

//	.....

}

int main(int argc, char *argv[])
{
	FILE	*samplefile;
	char	sampledirname[MAXLEN] = "c:/tmp/sample0";
	char	samplefilename[MAXLEN];
	char	outputfilename[MAXLEN] = "c:/tmp/markov_output.txt";
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
		samplecount++;
		fclose(samplefile);
	}

	// Initial setting
	printf("Init a , b and pi\n");
	for (i = 0; i < MAXSTATE; i++)
	{
		for (j = 0; j < MAXSTATE; j++)
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
		for (i = 0; i < MAXSTATE; i++)
		{
			// update a
			for (j = 0; j < MAXSTATE; j++)
			{
				na[i][j] = 0;
				bunbo = 0;
				for (n = 0; n < MAXSAMPLE; n++)
				{
					for (k = 0; k < maxlen[n]; k++)
					{
						na[i][j] += c[n][k][i][j];
						for (l = 0; l < MAXSTATE; l++)
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
				for (j = 0; j < MAXSTATE; j++)
					npi[i] += c[n][0][i][j];
			}

			// update b
			for (j = 0; j < MAXSYMBOL; j++)
			{
				nb[i][j] = 0;
				bunbo = 0;
				for (n = 0; n < MAXSAMPLE; n++)
				{
					for (k = 0; k < maxlen[n]; k++)
					{
						if (o[n][k] == j)
							for (l = 0; l < MAXSTATE; l++)
								nb[i][j] += c[n][k][i][l];
						for (l = 0; l < MAXSTATE; l++)
							bunbo += c[n][k][i][l];
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
		for (bunbo = 0, i = 0; i < MAXSTATE; i++)
		{
			bunbo += npi[i];
		}
		for (i = 0; i < MAXSTATE; i++)
			npi[i] /= bunbo;
		//	printf("continue -> ");
		for (i = 0; i < MAXSTATE; i++)
		{
			for (j = 0; j < MAXSTATE; j++)
			{
				err = fabs(a[i][j] - na[i][j]);
				if (err / a[i][j] > ERROR)
				{
					printf("error on a[%d][%d] = %g %g > %g\n", i, j, err, err / a[i][j], ERROR);
					break;
				}
			}
			if (j < MAXSTATE)
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
		if (!(i < MAXSTATE))
		{
			converg = 1;
		}
		else
		{
			converg = 0;
		}
		for (i = 0; i < MAXSTATE; i++)
		{
			for (j = 0; j < MAXSTATE; j++)
			{
				a[i][j] = na[i][j];
			}
			for (j = 0; j < MAXSYMBOL; j++)
			{
				b[i][j] = nb[i][j];
			}
			pi[i] = npi[i];
		}
	}
	printf("loop end and print HMM\n");
	HMMprint(outputfilename);
}
