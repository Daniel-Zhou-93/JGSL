#include <jni.h>
#include <assert.h>
#include <stdlib.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

JNIEXPORT void JNICALL Java_JGSL_SingularValueDecomposition_dcomp// A=U*S*(V)^T
  (JNIEnv *env, jclass obj, jint m, jint n, jint l, jdoubleArray ja, 
  jdoubleArray jv, jdoubleArray js){
      
      double *aElems = (*env)-> GetDoubleArrayElements (env, ja, NULL);
      double *vElems = (*env)-> GetDoubleArrayElements (env, jv, NULL);
      double *sElems = (*env)-> GetDoubleArrayElements (env, js, NULL);

      gsl_matrix_view tempa = gsl_matrix_view_array(aElems, m, n);//output is U
      gsl_matrix_view tempv = gsl_matrix_view_array(vElems, n, n);//v is untransposed form
      gsl_vector_view temps = gsl_vector_view_array(sElems, l);//singular value
      gsl_vector *work = gsl_vector_alloc(l); 

      gsl_linalg_SV_decomp( &tempa.matrix, &tempv.matrix, &temps.vector, work);
      gsl_vector_free(work);
      
      (*env)-> ReleaseDoubleArrayElements (env, ja, aElems, 0);
      (*env)-> ReleaseDoubleArrayElements (env, js, sElems, 0);
      (*env)-> ReleaseDoubleArrayElements (env, jv, vElems, 0);
  }


JNIEXPORT void JNICALL Java_JGSL_SingularValueDecomposition_slve
  (JNIEnv *env, jclass obj, jint m, jint n, jint l, jdoubleArray ja,
   jdoubleArray jv, jdoubleArray js, jdoubleArray jb, jdoubleArray jx)
{

      double *aElems = (*env)-> GetDoubleArrayElements (env, ja, NULL);
      double *vElems = (*env)-> GetDoubleArrayElements (env, jv, NULL);
      double *sElems = (*env)-> GetDoubleArrayElements (env, js, NULL);
      double *bElems = (*env)-> GetDoubleArrayElements (env, jb, NULL);
      double *xElems = (*env)-> GetDoubleArrayElements (env, jx, NULL);

      gsl_matrix_view tempa = gsl_matrix_view_array(aElems, m, n);//output is U
      gsl_matrix_view tempv = gsl_matrix_view_array(vElems, n, n);//v is untransposed form
      gsl_vector_view temps = gsl_vector_view_array(sElems, l);//singular value
      gsl_vector_view tempb = gsl_vector_view_array(bElems, m);//right hand
      gsl_vector_view tempx = gsl_vector_view_array(xElems, n);//solution
      

      gsl_linalg_SV_solve(&tempa.matrix, &tempv.matrix, &temps.vector, &tempb.vector, &tempx.vector);
      
      (*env)-> ReleaseDoubleArrayElements (env, ja, aElems, 0);
      (*env)-> ReleaseDoubleArrayElements (env, js, sElems, 0);
      (*env)-> ReleaseDoubleArrayElements (env, jv, vElems, 0);
      (*env)-> ReleaseDoubleArrayElements (env, jb, bElems, 0);
      (*env)-> ReleaseDoubleArrayElements (env, jx, xElems, 0);
}

