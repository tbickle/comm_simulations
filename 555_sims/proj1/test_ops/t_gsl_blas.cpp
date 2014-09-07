#include <iostream>
#include <vector>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>
using namespace std;

int main (void)
{
// init
	unsigned int NUM = 10;
	gsl_vector_complex* cvec = gsl_vector_complex_calloc(NUM);
	gsl_complex result;
	gsl_complex c;


// complex vector ops
	// complex
	c = gsl_complex_rect(1.0,2.0);
	gsl_vector_complex_set_all(cvec,c);
	// blas
	gsl_blas_zdscal(1.5,cvec);
	gsl_blas_zdotu(cvec,cvec,&result);


// print
	for (int i=0; i<NUM; i++) {cout<< GSL_REAL(gsl_vector_complex_get(cvec,i)) << "+j" << GSL_IMAG(gsl_vector_complex_get(cvec,i)) << endl;}
	cout << "result = " << GSL_REAL(result) << "+i" << GSL_IMAG(result) << endl; // zdotu(1+2i) = -30+i40
	// complex
	//c = gsl_complex_add(c,c);
	//cout << "c = " << GSL_REAL(c) << "+i" << GSL_IMAG(c) << endl;


// end
	gsl_vector_complex_free(cvec);
	return 0;
}

