#include <stdio.h>
#include <iostream>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
using namespace std;

int main (void)
{
// init
	unsigned int NUM = 5;
	gsl_vector* vec = gsl_vector_calloc(NUM); // create & initialize vector w/ 0's
	gsl_vector* a2 = gsl_vector_calloc(NUM);
	gsl_complex c;

// calc
	// vector ops
	gsl_vector_add_constant(vec,1);
	gsl_vector_set(vec,0,1.5);
	gsl_vector_add_constant(a2,2);
	gsl_vector_scale(vec,3);
	gsl_vector_mul(vec,a2);
	// complex
	c = gsl_complex_rect(1.0,2.0);
	c = gsl_complex_add(c,c);


// print
	for (int i=0; i<NUM; i++) {cout<< gsl_vector_get(vec,i) << endl;}
	cout << GSL_REAL(c) << "+j" << GSL_IMAG(c) << endl;


// end
	gsl_vector_free(vec);
	gsl_vector_free(a2);
	return 0;
}

