#include <jni.h>
#include <assert.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_matrix.h>

JNIEXPORT void JNICALL Java_JGSL_EigenvalueDecomposition_nonsymmv
  (JNIEnv *env, jclass obj, jint n, jdoubleArray ja, jdoubleArray jreal, 
   jdoubleArray jimag, jdoubleArray jmreal, jdoubleArray jmimag)
{
      
    double *aElems = (*env)-> GetDoubleArrayElements (env, ja, NULL);
    double *real = (*env)-> GetDoubleArrayElements (env, jreal, NULL); //eigenvalue real
    double *imag = (*env)-> GetDoubleArrayElements (env, jimag, NULL); //eigenvalue image
    double *mreal = (*env)-> GetDoubleArrayElements (env, jmreal, NULL); //eigenvector real
    double *mimag = (*env)-> GetDoubleArrayElements (env, jmimag, NULL); //eigenvector complex
    
    gsl_eigen_nonsymmv_workspace *w = gsl_eigen_nonsymmv_alloc(n);
    gsl_matrix_view tempa = gsl_matrix_view_array(aElems, n, n);
    gsl_vector_complex *eval = gsl_vector_complex_alloc (n);
    gsl_matrix_complex *evec = gsl_matrix_complex_alloc (n, n);

    gsl_eigen_nonsymmv( &tempa.matrix, eval, evec, w);
    gsl_eigen_nonsymmv_free(w);

    gsl_eigen_nonsymmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);
    int i, j;
    for (i=0; i<n; i++){
        gsl_complex eval_i = gsl_vector_complex_get(eval, i);
        *(real+i) = eval_i.dat[0];
        *(imag+i) = eval_i.dat[1];
        
        gsl_vector_complex_view eveci = gsl_matrix_complex_row (evec, i);// row
        for (j=0; j<n; j++){
            gsl_complex evec_i = gsl_vector_complex_get( &eveci.vector, j);
            *(mreal+i*n+j) = evec_i.dat[0];
            *(mimag+i*n+j) = evec_i.dat[1];
        }
    }

    gsl_vector_complex_free(eval);
    gsl_matrix_complex_free(evec);

    (*env)-> ReleaseDoubleArrayElements (env, ja, aElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, jreal, real, 0);
    (*env)-> ReleaseDoubleArrayElements (env, jimag, imag, 0);
    (*env)-> ReleaseDoubleArrayElements (env, jmreal, mreal, 0);
    (*env)-> ReleaseDoubleArrayElements (env, jmimag, mimag, 0);
}

JNIEXPORT void JNICALL Java_JGSL_EigenvalueDecomposition_symmv
  (JNIEnv *env, jclass obj, jint n, jdoubleArray ja,  jdoubleArray jeval, jdoubleArray jevec)
{
    double *aElems = (*env)-> GetDoubleArrayElements (env, ja, NULL);
    double *eval = (*env)-> GetDoubleArrayElements (env, jeval, NULL);
    double *evec = (*env)-> GetDoubleArrayElements (env, jevec, NULL);
    
    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(n);
    gsl_matrix_view tempa = gsl_matrix_view_array(aElems, n, n);
    gsl_vector_view tempeval = gsl_vector_view_array(eval, n); //eigenvector unordered
    gsl_matrix_view tempevec = gsl_matrix_view_array(evec, n, n); //corresponding eigenvectors

    gsl_eigen_symmv(&tempa.matrix, &tempeval.vector, &tempevec.matrix, w);
    gsl_eigen_symmv_free(w);
    
    gsl_eigen_symmv_sort( &tempeval.vector, &tempevec.matrix, GSL_EIGEN_SORT_ABS_ASC);//ascending order in magnitude

    (*env)-> ReleaseDoubleArrayElements (env, ja, aElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, jeval, eval, 0);
    (*env)-> ReleaseDoubleArrayElements (env, jevec, evec, 0);
}

