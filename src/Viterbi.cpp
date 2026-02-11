#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define	MAXSAMPLE	1000
#define	MAXSTATE	6    // Change 4 => 6.
#define	MAXSYMBOL	10
#define	MAXLEN		512

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
int	maxlen[MAXSAMPLE];		/* �T���v���t�@�C�����Ƃ̒������i�[ */
double	a[MAXSTATE][MAXSTATE];
double	b[MAXSTATE][MAXSYMBOL];
double	pi[MAXSTATE];
double	delta[MAXSAMPLE][MAXLEN][MAXSTATE];
int	phi[MAXSAMPLE][MAXLEN][MAXSTATE];
int	o[MAXSAMPLE][MAXLEN];
int	state[MAXSAMPLE][MAXLEN];
int states;

int HMMload(char *mn)
{
	int	i, j, f, k;
	FILE	*markovfile;
	char	lognum[MAXLEN];
	char	markovname[MAXLEN];
	char	c;
	float	loadval;
	char	tmp[MAXLEN];

	printf("hmmload start\n");
	strncpy(markovname, mn, MAXLEN);
	if ((markovfile = fopen(markovname, "r")) == NULL)
	{
		printf("markov file read open error : %s\n", markovname);
		exit(0);
	}
	for (i = 0; i < MAXSTATE; i++)
	{
		k = fscanf(markovfile, "%s", tmp);
		if (k == EOF && i > 0)
			break;
		if (k != 1 || strncmp(tmp, "State", 5))
		{
			printf("hmm file format error : State <-> %s\n", tmp);
			exit(0);
		}
		if (fscanf(markovfile, "%s", tmp) != 1)
		{
			printf("hmm file format error : tmp = %s\n", tmp);
			exit(0);
		}
		f = tmp[0] - 'A';
		if (fscanf(markovfile, "%g", &loadval) != 1)
		{
			printf("hmm file format error : loadval = %g\n", loadval);
			exit(0);
		}
		pi[f] = loadval;
		printf("f = %d , pi[%d] = %g\n", f, f, pi[f]);
		if (fscanf(markovfile, "%s", tmp) != 1 || strncmp(tmp, "Output", 6))
		{
			printf("hmm file format error : Output <-> %s\n", tmp);
			exit(0);
		}
		printf("Output: loading\n");
		for (j = 0; j < MAXSYMBOL; j++)
		{
			if (fscanf(markovfile, "%g", &loadval) != 1)
				break;
			b[f][j] = loadval;
			// printf("b[%d][%d] = %g\n", f, j, b[f][j]);
		}
		if (fscanf(markovfile, "%s", tmp) != 1 && strncmp(tmp, "Transition", 10))
		{
			printf("hmm file format error : Transition <-> %s\n", tmp);
			exit(0);
		}
		printf("Transition: loading\n");
		for (j = 0; j < MAXSTATE; j++)
		{
			if (fscanf(markovfile, "%g", &loadval) != 1)
				break;
			a[f][j] = loadval;
			// printf("a[%d][%d] = %g\n", f, j, a[f][j]);
		}
	}
	printf("load loop end with %d State\n", i);
	fclose(markovfile);

	printf("loaded HMM is\n");
	for (i = 0; i < MAXSTATE; i++)
	{
		printf("State %c %g\n", i + 'A', pi[i]);
		printf("Output");
		for (j = 0; j < MAXSYMBOL; j++)
		{
			printf(" %g", b[i][j]);
		}
		printf("\n");
		printf("Transition");
		for (j = 0; j < MAXSTATE; j++)
		{
			printf(" %g", a[i][j]);
		}
		printf("\n");
	}
	return 0;
}

double viterbi(int samplenum)
{
	int	i, j, k, l;
	double	maxval;
	int		maxindex;
	double	outputprob;

	for (i = 0; i < MAXSTATE; i++)
	{
		delta[samplenum][0][i] = pi[i] * b[i][o[samplenum][0]];
		phi[samplenum][0][i] = 0;
	}
	for (i = 1; i < maxlen[samplenum]; i++)
	{
		for (j = 0; j < MAXSTATE; j++)
		{
//			.....
			maxval = -1.0;
            maxindex = 0;
            for (k = 0; k < states; k++) {
                // Calculate the probability of transitioning from the previous state k to the current state j.
                double val = delta[samplenum][i - 1][k] * a[k][j];
                if (val > maxval) {
                    maxval = val;
                    maxindex = k;
                }
            }
            // Update delta with the maximum path probability and the current output probability.
            delta[samplenum][i][j] = maxval * b[j][o[samplenum][i]];
            // Store the index of the best previous state to allow for backtracking.
            phi[samplenum][i][j] = maxindex;
		}
	}
	for (outputprob = 0, k = 0; k < states; k++)
	{
		if (delta[samplenum][maxlen[samplenum] - 1][k] > outputprob)
		{
			outputprob = delta[samplenum][maxlen[samplenum] - 1][k];
			state[samplenum][maxlen[samplenum] - 1] = k;
		}
	}
	for (k = maxlen[samplenum] - 1; k > 0; k--)
		state[samplenum][k - 1] = phi[samplenum][k][state[samplenum][k]];

	return(outputprob);
}

int main(int argc, char *argv[])
{
	// input model file and others.
	if (argc < 6) {
    	printf("Usage: %s [SampleName] [OutputDir] [TargetModelFile] [States] [Seed]\n", argv[0]);
    	exit(1);
	}
	char *sampledirname = argv[1];
	char *outputdir = argv[2];
	char *markovfilename = argv[3];
	states = atoi(argv[4]);
	int seed = atoi(argv[5]);
 	// Use these variables to construct file paths later.

	FILE	*samplefile;
//	char	sampledirname[MAXLEN] = "c:/tmp/sample0";
	char	samplefilename[MAXLEN];
//	char	markovfilename[MAXLEN] = "c:/tmp/markov_output0.txt";
	char	tmp[MAXLEN];

	int	i, j, k, l, m, n;
	int	readsymbol;
	int	samplecount, timecount;
	double	outputprob;

	for (samplecount = 0, i = 0; i < MAXSAMPLE; i++)
	{
	//	strncpy(samplefilename, sampledirname, MAXLEN);
	//	strncat(samplefilename, "/", MAXLEN);
	//	sprintf(tmp, "%d", i + 1);
	//	strncat(samplefilename, tmp, MAXLEN);
	//	strncat(samplefilename, ".txt", MAXLEN);
    	sprintf(samplefilename, "%s/%d.txt", sampledirname, i + 1);    
		// The path format: [OutputDir]/[Index].txt.

		if ((samplefile = fopen(samplefilename, "r")) == NULL)
			continue;
		printf("file open success : %s\n", samplefilename);
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
	HMMload(markovfilename);


	for (i = 0; i < samplecount; i++)
	{
		outputprob = viterbi(i);
		// strncpy(samplefilename, sampledirname, MAXLEN);
		// strncat(samplefilename, "/", MAXLEN);
		// sprintf(tmp, "%d", i + 1)
		// strncat(samplefilename, tmp, MAXLEN);
		// strncat(samplefilename, "_out.txt", MAXLEN);	// samplefilename = "(sampledirname)/1234_out.txt"
		sprintf(samplefilename, "%s/s%d_%d_out.txt", outputdir, states, i + 1); // The path format: [OutputDir]/s[States]_[Index]_out.txt.

		if ((samplefile = fopen(samplefilename, "w")) == NULL)
			continue;
		for (j = 0; j < maxlen[i]; j++)
		{
			if (j > 0)
				fprintf(samplefile, " ");
			fprintf(samplefile, "%c", 'A' + state[i][j]);
		}
		fclose(samplefile);
		printf("outputprob = %g , file: %s\n", outputprob, samplefilename);
	}
}
