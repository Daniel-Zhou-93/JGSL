#include <jni.h>
#include <assert.h>
#include <stdlib.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

JNIEXPORT int JNICALL Java_JGSL_CholeskyDecomposition_dcomp
  (JNIEnv *env, jclass obj, jint n, jdoubleArray ja, jint info)
{
    double *aElems = (*env)-> GetDoubleArrayElements (env, ja, NULL);
    gsl_matrix_view tempa=gsl_matrix_view_array(aElems, n, n);
    info=gsl_linalg_cholesky_decomp1(&tempa.matrix);
    
    (*env)-> ReleaseDoubleArrayElements (env, ja, aElems, 0);
    return info;
}


JNIEXPORT void JNICALL Java_JGSL_CholeskyDecomposition_invert
  (JNIEnv *env, jclass obj, jint n, jdoubleArray ja)
{
    double *aElems = (*env)-> GetDoubleArrayElements (env, ja, NULL);

    gsl_matrix_view tempa=gsl_matrix_view_array(aElems, n, n);
    gsl_linalg_cholesky_invert(&tempa.matrix);

    (*env)-> ReleaseDoubleArrayElements (env, ja, aElems, 0);
}


JNIEXPORT void JNICALL Java_JGSL_CholeskyDecomposition_slve
  (JNIEnv *env, jclass obj, jint n, jdoubleArray ja, jdoubleArray jb, jdoubleArray jx)
{
    double *aElems = (*env)-> GetDoubleArrayElements (env, ja, NULL);
    double *bElems = (*env)-> GetDoubleArrayElements (env, jb, NULL);
    double *xElems = (*env)-> GetDoubleArrayElements (env, jx, NULL);

    gsl_vector_view tempb=gsl_vector_view_array(bElems, n);
    gsl_vector_view tempx=gsl_vector_view_array(xElems, n);
    gsl_matrix_view tempa= gsl_matrix_view_array(aElems, n, n);

    gsl_linalg_cholesky_solve(&tempa.matrix, &tempb.vector, &tempx.vector);

    (*env)-> ReleaseDoubleArrayElements (env, ja, aElems, JNI_ABORT);
    (*env)-> ReleaseDoubleArrayElements (env, jb, bElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, jx, xElems, 0);
}
