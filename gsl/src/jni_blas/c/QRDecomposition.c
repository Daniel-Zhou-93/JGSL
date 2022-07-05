#include <jni.h>
#include <assert.h>
#include <stdlib.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

JNIEXPORT void JNICALL Java_JGSL_QRDecomposition_dcomp
  (JNIEnv *env, jclass obj, jint m, jint n, jint l, jdoubleArray ja, jdoubleArray jtau)
{
    double *aElems = (*env)-> GetDoubleArrayElements (env, ja, NULL);
    double *tauElems = (*env)-> GetDoubleArrayElements (env, jtau, NULL);
    gsl_vector_view temptau=gsl_vector_view_array(tauElems, l);
    gsl_matrix_view tempa= gsl_matrix_view_array(aElems, m, n);

    gsl_linalg_QR_decomp(&tempa.matrix, &temptau.vector);

      
    (*env)-> ReleaseDoubleArrayElements (env, ja, aElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, jtau, tauElems, 0);
}


JNIEXPORT void JNICALL Java_JGSL_QRDecomposition_unpack
  (JNIEnv *env, jclass obj, jint m, jint n, jint l,
   jdoubleArray ja, jdoubleArray jtau, jdoubleArray jq)
{
    double *aElems = (*env)-> GetDoubleArrayElements (env, ja, NULL);
    double *tauElems = (*env)-> GetDoubleArrayElements (env, jtau, NULL);
    double *qElems = (*env)-> GetDoubleArrayElements (env, jq, NULL);
    gsl_matrix *r= gsl_matrix_alloc( m, n);

    gsl_vector_view temptau = gsl_vector_view_array(tauElems, l);
    gsl_matrix_view tempa = gsl_matrix_view_array(aElems, m, n);
    gsl_matrix_view tempq = gsl_matrix_view_array(qElems, m, m);

    gsl_linalg_QR_unpack( &tempa.matrix, &temptau.vector, &tempq.matrix, r);
    gsl_matrix_free(r);
    (*env)-> ReleaseDoubleArrayElements (env, ja, aElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, jtau, tauElems, JNI_ABORT);
    (*env)-> ReleaseDoubleArrayElements (env, jq, qElems, 0);
}


JNIEXPORT void Java_JGSL_QRDecomposition_slve
(JNIEnv *env, jclass klass, jint m, jint n, jint l, jdoubleArray ja, jdoubleArray jtau
 , jdoubleArray jb, jdoubleArray jx){
    
    double *aElems = (*env)-> GetDoubleArrayElements (env, ja, NULL);
    double *bElems = (*env)-> GetDoubleArrayElements (env, jb, NULL);
    double *tauElems = (*env)-> GetDoubleArrayElements (env, jtau, NULL);
    double *xElems = (*env)-> GetDoubleArrayElements (env, jx, NULL);    
    assert(aElems && tauElems && bElems && xElems);

    gsl_vector_view temptau = gsl_vector_view_array(tauElems, l);
    gsl_vector_view tempb = gsl_vector_view_array(bElems, m);
    gsl_matrix_view tempa = gsl_matrix_view_array(aElems, m, n);
    gsl_vector_view tempx = gsl_vector_view_array(xElems, n);

    gsl_linalg_QR_solve( &tempa.matrix, &temptau.vector, &tempb.vector, &tempx.vector); 

    
    (*env)-> ReleaseDoubleArrayElements (env, ja, aElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, jtau, tauElems, JNI_ABORT);
    (*env)-> ReleaseDoubleArrayElements (env, jb, bElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, jx, xElems, 0);
}


