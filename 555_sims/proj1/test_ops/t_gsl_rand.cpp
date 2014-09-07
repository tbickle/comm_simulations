#include <stdio.h>
#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>

/*
int main (void)
{
	const gsl_rng_type* T;
	gsl_rng* r;

	int i, n = 10;
	//double mu = 3.0;
	double sigma = 1.0;

	// create a generator chosen by the environment variable GSL_RNG_TYPE
	gsl_rng_env_setup();

	T = gsl_rng_default;
	r = gsl_rng_alloc (T);

	// print n random variates chosen from the poisson distribution with mean parameter mu
	for (i=0; i<n; i++) {
		//unsigned int k = gsl_ran_poisson (r, mu);
		double k = gsl_ran_gaussian (r, sigma);
		printf (" %f", k);
	}

	printf ("\n");
	gsl_rng_free (r);
	return 0;
}
*/

gsl_rng* r;  // global generator

int
main (void)
{
  const gsl_rng_type* T;

  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  
  printf ("generator type: %s\n", gsl_rng_name (r));
  printf ("seed = %lu\n", gsl_rng_default_seed);
  printf ("first value = %lu\n", gsl_rng_get (r));

  gsl_rng_free (r);
  return 0;
}
