#include <jni.h>
#include <assert.h>
#include <stdlib.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>


/* LU */
JNIEXPORT int JNICALL Java_JGSL_LUDecomposition_decomp
  (JNIEnv *env, jclass obj, jint n, jdoubleArray ja, jlongArray jp, jint signum)
{
    double *aElems = (*env)-> GetDoubleArrayElements (env, ja, NULL);
    long *pElems = (*env)-> GetLongArrayElements (env, jp, NULL);

    gsl_permutation *p= gsl_permutation_alloc(n);
    gsl_matrix_view tempa=gsl_matrix_view_array(aElems, n, n);

    gsl_linalg_LU_decomp(&tempa.matrix, p, &signum);
    
    //int i;
    //for(i=0; i<n ; i++)
    //    *(pElems+i)=(long) *(p->data+i);

    pElems= (long *) p->data;
    //gsl_permutation_free (p);
    (*env)-> ReleaseDoubleArrayElements (env, ja, aElems, 0);
    (*env)-> ReleaseLongArrayElements (env, jp, pElems, 0);

    return signum;
}

JNIEXPORT void JNICALL Java_JGSL_LUDecomposition_slve
  (JNIEnv *env, jclass obj, jint n, jdoubleArray ja, jdoubleArray jb, 
   jlongArray jp, jdoubleArray jx)
{
    double *aElems = (*env)-> GetDoubleArrayElements (env, ja, NULL);
    double *bElems = (*env)-> GetDoubleArrayElements (env, jb, NULL);
    double *xElems = (*env)-> GetDoubleArrayElements (env, jx, NULL);
    long *pElems = (*env)-> GetLongArrayElements (env, jp, NULL);

    gsl_permutation *p= gsl_permutation_alloc(n);
    p->data=(size_t *) pElems;
    gsl_vector_view tempb=gsl_vector_view_array(bElems, n);
    gsl_vector_view tempx=gsl_vector_view_array(xElems, n);
    gsl_matrix_view tempa= gsl_matrix_view_array(aElems, n, n);

    gsl_linalg_LU_solve(&tempa.matrix, p, &tempb.vector, &tempx.vector);

    //int i;
    //for(i=0; i<n ; i++)
    //    *(pElems+i)=(long) *(p->data+i);
    //gsl_permutation_free (p);
    pElems= (long *) p->data;


    (*env)-> ReleaseDoubleArrayElements (env, ja, aElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, jb, bElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, jx, xElems, 0);
    (*env)-> ReleaseLongArrayElements (env, jp, pElems, 0);
    //gsl_permutation_free(p);
}

JNIEXPORT void JNICALL Java_JGSL_LUDecomposition_invert
  (JNIEnv *env, jclass obj, jint n, jdoubleArray ja, jdoubleArray jb, 
   jlongArray jp)
{
    double *aElems = (*env)-> GetDoubleArrayElements (env, ja, NULL);
    double *bElems = (*env)-> GetDoubleArrayElements (env, jb, NULL);
    long *pElems = (*env)-> GetLongArrayElements (env, jp, NULL);

    gsl_permutation *p= gsl_permutation_alloc(n);
    p->data=(size_t *) pElems;
    gsl_matrix_view tempb=gsl_matrix_view_array(bElems, n, n);
    gsl_matrix_view tempa= gsl_matrix_view_array(aElems, n, n);

    gsl_linalg_LU_invert (&tempa.matrix, p, &tempb.matrix);

    //int i;
    //for(i=0; i<n ; i++)
    //    *(pElems+i)=(long) *(p->data+i);
    //gsl_permutation_free (p);
    pElems= (long *) p->data;


    (*env)-> ReleaseDoubleArrayElements (env, ja, aElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, jb, bElems, 0);
    (*env)-> ReleaseLongArrayElements (env, jp, pElems, 0);
    //gsl_permutation_free(p);
}
