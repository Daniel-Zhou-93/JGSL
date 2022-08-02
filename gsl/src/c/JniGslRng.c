#define MATHLIB_STANDALONE

#include <jni.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
//univariate multivariate, normal, gamma, inverse gamma, wishart, inverse wishart
        //extern const gsl_rng_type *gsl_rng_rand;
	//extern unsigned long int gsl_rng_default_seed;
	//T = gsl_rng_rand;
        gsl_rng *r;

        //const gsl_rng_type *T = gsl_rng_rand;
	//T = gsl_rng_rand;
        //r = gsl_rng_alloc (T);
JNIEXPORT void JNICALL Java_JGSL_JniGslRng_seed
    (JNIEnv *env, jclass obj)
{

	const gsl_rng_type *T;
	T = gsl_rng_rand;
	r = gsl_rng_alloc(T);
	//gsl_rng_set(r, time(NULL));
        //gsl_rng_env_setup();
}

JNIEXPORT jdouble JNICALL Java_JGSL_JniGslRng_uninorm
    (JNIEnv *env, jclass obj, jdouble sigma)
{

	//const gsl_rng_type *T;
	//T = gsl_rng_rand;
	//r = gsl_rng_alloc(T);
	//gsl_rng_set(r, time(NULL));
        gsl_rng_env_setup();
        double result;
	result = gsl_ran_gaussian(r, sigma);
        //gsl_rng_free (r);
        return result;
}

JNIEXPORT void JNICALL Java_JGSL_JniGslRng_multinorm
  (JNIEnv *env, jclass obj, jint k, jdoubleArray jmu, jdoubleArray jl, jdoubleArray jresult)
{
    double *muElems = (*env)-> GetDoubleArrayElements (env, jmu, NULL);
    double *lElems = (*env)-> GetDoubleArrayElements (env, jl, NULL);
    double *resultElems = (*env)-> GetDoubleArrayElements (env, jresult, NULL);
    //const gsl_rng_type *T;
    //gsl_rng *r;
    //T = gsl_rng_rand;
    //r = gsl_rng_alloc (T);
    //gsl_rng_set(r, time(NULL));
    gsl_rng_env_setup();
    gsl_vector_view tempmu = gsl_vector_view_array(muElems, k);
    gsl_matrix_view templ = gsl_matrix_view_array(lElems, k, k);
    gsl_linalg_cholesky_decomp(&templ.matrix);
    gsl_vector_view tempresult = gsl_vector_view_array(resultElems, k);

    gsl_ran_multivariate_gaussian( r, &tempmu.vector, &templ.matrix, &tempresult.vector);
    //gsl_rng_free (r);
    (*env)-> ReleaseDoubleArrayElements (env, jmu, muElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, jl, lElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, jresult, resultElems, 0);
}

JNIEXPORT jdouble JNICALL Java_JGSL_JniGslRng_gamma
    (JNIEnv *env, jclass obj, jdouble a, jdouble b)
{
        //const gsl_rng_type *T;
        //gsl_rng *r;
        //T = gsl_rng_rand;
        //r = gsl_rng_alloc (T);
	//gsl_rng_set(r, time(NULL));
        gsl_rng_env_setup();
        double result;

	result = gsl_ran_gamma(r, a, b);
        //gsl_rng_free (r);
        return result;
}

JNIEXPORT jdouble JNICALL Java_JGSL_JniGslRng_invergamma
    (JNIEnv *env, jclass obj, jdouble a, jdouble b)
{
        //const gsl_rng_type *T;
        //gsl_rng *r;
        //T = gsl_rng_rand;
        //r = gsl_rng_alloc (T);
	//gsl_rng_set(r, time(NULL));
        gsl_rng_env_setup();
        double result;

	result = 1/gsl_ran_gamma(r, a, 1/b);
        //gsl_rng_free (r);
        return result;
}

//NOTE: Pass in the lower triangular of the scale matrix for the Wishart, and the lower triangular for the inverse of the matrix
//      for the inverse Wishart to achieve the desired behavior.

JNIEXPORT void JNICALL Java_JGSL_JniGslRng_wishart
  (JNIEnv *env, jclass obj, jdouble n, jint p, jdoubleArray jl, jdoubleArray jresult)
{
    double *lElems = (*env)-> GetDoubleArrayElements (env, jl, NULL);
    double *resultElems = (*env)-> GetDoubleArrayElements (env, jresult, NULL);
    //const gsl_rng_type *T;
    gsl_matrix *work = gsl_matrix_alloc(p, p);
    //gsl_rng *r;
    //T = gsl_rng_rand;
    //r = gsl_rng_alloc (T);
    //gsl_rng_set(r, time(NULL));
    gsl_rng_env_setup();
    gsl_matrix_view templ = gsl_matrix_view_array(lElems, p, p);
    gsl_matrix_view tempresult = gsl_matrix_view_array(resultElems, p, p);

    gsl_ran_wishart( r, n, &templ.matrix, &tempresult.matrix, work);
    //gsl_rng_free (r);
    gsl_matrix_free(work);

    (*env)-> ReleaseDoubleArrayElements (env, jl, lElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, jresult, resultElems, 0);
}


JNIEXPORT void JNICALL Java_JGSL_JniGslRng_inverwishart
  (JNIEnv *env, jclass obj, jdouble n, jint p, jdoubleArray jl, jdoubleArray jresult)
{
    double *lElems = (*env)-> GetDoubleArrayElements (env, jl, NULL);
    double *resultElems = (*env)-> GetDoubleArrayElements (env, jresult, NULL);
    //const gsl_rng_type *T;
    int signum=1;
    gsl_matrix *work = gsl_matrix_alloc(p, p);
    gsl_matrix *inverse = gsl_matrix_alloc(p, p);
    //gsl_rng *r;

    //T = gsl_rng_rand;
    //r = gsl_rng_alloc (T);
    //gsl_rng_set(r, time(NULL));
    gsl_rng_env_setup();
    gsl_matrix_view templ = gsl_matrix_view_array(lElems, p, p);
    gsl_matrix_view tempresult = gsl_matrix_view_array(resultElems, p, p);

    gsl_ran_wishart( r, n, &templ.matrix, inverse, work);
    //gsl_rng_free (r);
    gsl_matrix_free(work);

    gsl_permutation *q= gsl_permutation_alloc(p);

    gsl_linalg_LU_decomp(inverse, q, &signum);
    gsl_linalg_LU_invert (inverse, q, &tempresult.matrix);
    gsl_permutation_free(q);
    gsl_matrix_free(inverse);

    (*env)-> ReleaseDoubleArrayElements (env, jl, lElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, jresult, resultElems, 0);
}

JNIEXPORT void JNICALL Java_JGSL_JniGslRng_free
    (JNIEnv *env, jclass obj)
{

	gsl_rng_free(r);
}

