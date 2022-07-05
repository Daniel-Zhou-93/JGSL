#include <jni.h>
#include <assert.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>

#define jniRowMajor 101
#define jniColMajor 102

#define jniNoTrans 111
#define jniTrans   112
#define jniConjTrans    113

#define jniUpper   121
#define jniLower   122

#define jniNonUnit 131
#define jniUnit    132

#define jniLeft    141
#define jniRight   142

/* Level 1: dscal, daxpy, ddot */

JNIEXPORT void Java_JGSL_Matrix_dscal
(JNIEnv *env, jclass klass, jint n, jdouble alpha, jdoubleArray x){
    
    /* dscal:  x = alpha * x */
    
    double *xElems;
    xElems = (*env)-> GetDoubleArrayElements (env, x, NULL);
    assert(xElems);
    
    gsl_vector_view tempx=gsl_vector_view_array(xElems,n);
    gsl_blas_dscal(alpha, &tempx.vector);

    (*env)-> ReleaseDoubleArrayElements (env, x, xElems, 0);
}


JNIEXPORT void Java_JGSL_Matrix_daxpy
(JNIEnv *env, jclass klass, jint n, jdouble alpha, jdoubleArray x,
 jdoubleArray y){
    
    /* daxpy: y = alpha * x + y */
    
    double *xElems, *yElems;
    xElems = (*env)-> GetDoubleArrayElements (env, x, NULL);
    yElems = (*env)-> GetDoubleArrayElements (env, y, NULL);
    assert(xElems && yElems);
    
    gsl_vector_view tempx=gsl_vector_view_array(xElems,n);
    gsl_vector_view tempy=gsl_vector_view_array(yElems,n);
    gsl_blas_daxpy(alpha, &tempx.vector, &tempy.vector);
    
    (*env)-> ReleaseDoubleArrayElements (env, y, yElems, 0); 
    (*env)-> ReleaseDoubleArrayElements (env, x, xElems, JNI_ABORT);
}

JNIEXPORT jdouble Java_JGSL_Matrix_ddot
(JNIEnv *env, jclass klass, jint n, jdoubleArray x, jdoubleArray y){
    
    /* ddot:  forms the dot product of two vectors x and y.*/
    
    double *xElems, *yElems;
    double result;
    
    xElems = (*env)-> GetDoubleArrayElements (env, x, NULL);
    yElems = (*env)-> GetDoubleArrayElements (env, y, NULL);
    assert(xElems && yElems);
    
    gsl_vector_view tempx=gsl_vector_view_array(xElems,n);
    gsl_vector_view tempy=gsl_vector_view_array(yElems,n);
    
    gsl_blas_ddot( &tempx.vector, &tempy.vector, &result);
    
    (*env)-> ReleaseDoubleArrayElements (env, y, yElems, JNI_ABORT);
    (*env)-> ReleaseDoubleArrayElements (env, x, xElems, JNI_ABORT);
    
    return result;
}

/* Level 2: dgemv, dtrmv, dsymv */

JNIEXPORT void Java_JGSL_Matrix_dgemv
(JNIEnv *env, jclass klass, jint Trans, jint m, jint n,
 jdouble alpha, jdoubleArray A, jdoubleArray x, jdouble beta, jdoubleArray y){
    
    /* DGEMV  performs one of the matrix-vector operations
     y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y, */
    
    double *AElems, *xElems, *yElems;
    
    AElems = (*env)-> GetDoubleArrayElements (env, A, NULL);
    xElems = (*env)-> GetDoubleArrayElements (env, x, NULL);
    yElems = (*env)-> GetDoubleArrayElements (env, y, NULL);
    assert(AElems && xElems && yElems);

    gsl_vector_view tempx=gsl_vector_view_array(xElems,n);
    gsl_vector_view tempy=gsl_vector_view_array(yElems,m);

    gsl_matrix_view tempA=gsl_matrix_view_array(AElems,m,n);
    if (Trans == jniNoTrans) {
        gsl_blas_dgemv(jniNoTrans, alpha, &tempA.matrix, &tempx.vector, beta, &tempy.vector);}
    else if (Trans == jniTrans) {
        gsl_blas_dgemv(jniTrans, alpha, &tempA.matrix, &tempx.vector, beta, &tempy.vector);}
    else if (Trans == jniConjTrans) {
        gsl_blas_dgemv(jniConjTrans, alpha, &tempA.matrix, &tempx.vector, beta, &tempy.vector);}
    else {fprintf(stderr, "** Illegal Trans setting \n"); return;}
    
    (*env)-> ReleaseDoubleArrayElements (env, y, yElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, x, xElems, JNI_ABORT);
    (*env)-> ReleaseDoubleArrayElements (env, A, AElems, JNI_ABORT);
}

JNIEXPORT void Java_JGSL_Matrix_dtrmv
(JNIEnv *env, jclass klass, jint Uplo, jint Trans, jint Diag,
 jint n, jdoubleArray A, jdoubleArray x){
    
    /*  DTRMV  performs one of the matrix-vector operations
     x := A*x,   or   x := A**T*x,
     where x is an n element vector and  A is an n by n unit, or non-unit,
     upper or lower triangular matrix. */
    
    // LDA specifies the first dimension of A; lda = n
    double *AElems, *xElems;
    int Ts, uplo, diag;
    
    AElems = (*env)-> GetDoubleArrayElements (env, A, NULL);
    xElems = (*env)-> GetDoubleArrayElements (env, x, NULL);    
    assert(AElems && xElems);    

    gsl_vector_view tempx=gsl_vector_view_array(xElems,n);

    
    if (Diag == jniNonUnit) {diag = jniNonUnit;}
    else if (Diag == jniUnit) {diag = jniUnit;}
    else {fprintf(stderr, "** Illegal Diag setting \n"); return;}

    gsl_matrix_view tempA=gsl_matrix_view_array(AElems,n,n);        
    if (Trans == jniNoTrans) {Ts = jniNoTrans;}
    else if (Trans == jniTrans) {Ts = jniTrans;}
    else if (Trans == jniConjTrans) {Ts = jniConjTrans;}
    else {fprintf(stderr, "** Illegal Trans setting \n"); return;}
        
    if (Uplo == jniUpper) {uplo = jniUpper;}                // A is an upper triangular matrix.
    else if (Uplo == jniLower) {uplo = jniLower;}           // A is a lower triangular matrix
    else {fprintf(stderr, "** Illegal Uplo setting \n"); return;}
    gsl_blas_dtrmv(uplo, Ts, diag, &tempA.matrix, &tempx.vector);        
    
    (*env)-> ReleaseDoubleArrayElements (env, x, xElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, A, AElems, JNI_ABORT);
}

JNIEXPORT void Java_JGSL_Matrix_dsymv
(JNIEnv *env, jclass klass, jint Uplo, jint n, jdouble alpha,
 jdoubleArray A, jdoubleArray x, jdouble beta, jdoubleArray y){
    
    /*  DSYMV  performs the matrix-vector  operation
     y := alpha*A*x + beta*y
     where alpha and beta are scalars, x and y are n element vectors and
     A is an n by n symmetric matrix */
    double *AElems, *xElems, *yElems;

    AElems = (*env)-> GetDoubleArrayElements (env, A, NULL);
    xElems = (*env)-> GetDoubleArrayElements (env, x, NULL);
    yElems = (*env)-> GetDoubleArrayElements (env, y, NULL);
    
    assert(AElems && xElems && yElems);

    gsl_vector_view tempx=gsl_vector_view_array(xElems,n);
    gsl_vector_view tempy=gsl_vector_view_array(yElems,n);
    gsl_matrix_view tempA=gsl_matrix_view_array(AElems,n,n);
        
    // LDA specifies the first dimension of A; lda = n

    if (Uplo == jniUpper) {
        gsl_blas_dsymv(CblasUpper, alpha, &tempA.matrix, &tempx.vector, beta, &tempy.vector);} // A is an upper triangular matrix.
    else if (Uplo == jniLower) {
        gsl_blas_dsymv(CblasLower, alpha,  &tempA.matrix, &tempx.vector, beta, &tempy.vector);} // A is a lower triangular matrix
    else {fprintf(stderr, "** Illegal Uplo setting \n"); return;}
    
    (*env)-> ReleaseDoubleArrayElements (env, y, yElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, x, xElems, JNI_ABORT);
    (*env)-> ReleaseDoubleArrayElements (env, A, AElems, JNI_ABORT);
    
}


/* Level 3: dgemm, dtrmm, dsymm, dtrsm */

JNIEXPORT void Java_JGSL_Matrix_dgemm
(JNIEnv *env, jclass klass, jint TransA, jint TransB, jint m,
 jint n, jint l, jint k, jdouble alpha, jdoubleArray  A, jdoubleArray B,
 jdouble beta, jdoubleArray C){
    /* DGEMM  performs one of the matrix-matrix operations
     C := alpha*op( A )*op( B ) + beta*C,
     where  op( X ) is one of op( X ) = X   or   op( X ) = X**T,
     alpha and beta are scalars, and A, B and C are matrices, with op( A )
     an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix. */

    double *aElems, *bElems, *cElems;
    int TA, TB, nrow, bcol;
    aElems = (*env)-> GetDoubleArrayElements (env,A,NULL);
    bElems = (*env)-> GetDoubleArrayElements (env,B,NULL);
    cElems = (*env)-> GetDoubleArrayElements (env,C,NULL);
    assert(aElems && bElems && cElems);

    gsl_matrix_view tempa=gsl_matrix_view_array(aElems,m,n);
    gsl_matrix_view tempb=gsl_matrix_view_array(bElems,l,k);


    if (TransA == jniNoTrans) { TA = jniNoTrans; nrow = m;}
    else if (TransA == jniTrans) {TA = jniTrans; nrow = n;}
    else if (TransA == jniConjTrans) {TA = jniConjTrans; nrow = n;}
    else {fprintf(stderr, "** Illegal TransA setting \n"); return;}

    if (TransB == jniNoTrans) { TB = jniNoTrans; bcol = k;}
    else if (TransB == jniTrans) { TB = jniTrans; bcol = l;}
    else if (TransB == jniConjTrans) {TB = jniConjTrans; bcol = l;}
    else {fprintf(stderr, "** Illegal TransB setting \n"); return;}

    gsl_matrix_view tempc=gsl_matrix_view_array(cElems,nrow,bcol);               
    gsl_blas_dgemm(TA, TB, alpha, &tempa.matrix, &tempb.matrix, beta, &tempc.matrix);

    (*env)-> ReleaseDoubleArrayElements (env, C, cElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, B, bElems, JNI_ABORT);
    (*env)-> ReleaseDoubleArrayElements (env, A, aElems, JNI_ABORT);
}


JNIEXPORT void Java_JGSL_Matrix_dtrmm
(JNIEnv *env, jclass klass, jint Side, jint Uplo, jint TransA,
 jint Diag, jint m, jint n, jdouble alpha, jdoubleArray  A, jdoubleArray B){
    
    /* dtrmm: performs one of the matrix-matrix operations
     B := alpha*op( A )*B,   or   B := alpha*B*op( A )
     where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
     non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
     op( A ) = A   or   op( A ) = A**T.*/
    
    double *aElems, *bElems;
    int side, uplo, TA, diag, i, j;
    aElems = (*env)-> GetDoubleArrayElements (env,A, NULL);
    bElems = (*env)-> GetDoubleArrayElements (env,B, NULL);        
    assert(aElems && bElems);

    gsl_matrix_view tempa=gsl_matrix_view_array(aElems,m,n);
    gsl_matrix_view tempb=gsl_matrix_view_array(bElems,m,n);

    if (Diag == jniNonUnit) {diag = jniNonUnit;}
    else if (Diag == jniUnit) {diag = jniUnit;}
    else {fprintf(stderr, "** Illegal Diag setting \n"); return;}

    if (TransA == jniNoTrans) {TA = jniNoTrans;}
    else if (TransA == jniTrans) {TA = jniTrans;}
    else if (TransA == jniConjTrans) {TA = jniConjTrans;}
    else {fprintf(stderr, "** Illegal TransA setting \n"); return;}
    
    if (Uplo == jniUpper){ uplo = jniUpper;}                // B := alpha*op( A )*B
    else if(Uplo == jniLower){ uplo = jniLower;}
    else{fprintf(stderr, "** Illegal Uplo setting \n"); return;}

    if (Side == jniLeft){
        side = jniLeft; 
        gsl_matrix_view tempa=gsl_matrix_view_array(aElems,m,m);
        gsl_blas_dtrmm( side, uplo, TA, diag, alpha, &tempa.matrix, &tempb.matrix);      
    }        // B := alpha*op( A )*B
    else if(Side == jniRight){
        side = jniRight;
        gsl_matrix_view tempa=gsl_matrix_view_array(aElems,n,n);
        gsl_blas_dtrmm( side, uplo, TA, diag, alpha, &tempa.matrix, &tempb.matrix);
    }    // B := alpha*B*op( A )
    else{fprintf(stderr, "** Illegal Side setting \n"); return;}

    (*env)-> ReleaseDoubleArrayElements (env, B, bElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, A, aElems, JNI_ABORT);
}

JNIEXPORT void Java_JGSL_Matrix_dsymm
(JNIEnv *env, jclass klass, jint Side, jint Uplo,
 jint m, jint n, jdouble alpha, jdoubleArray  A, jdoubleArray B,
 jdouble beta, jdoubleArray C){
    
    /*DSYMM  performs one of the matrix-matrix operations:
     C := alpha*A*B + beta*C, or C := alpha*B*A + beta*C,
     where A is a symmetric matrix and  B and C are  m by n matrices.*/
    
    double *aElems, *bElems, *cElems;
    int lda, i, j;
    int side, uplo;
    aElems = (*env)-> GetDoubleArrayElements (env,A, NULL);
    bElems = (*env)-> GetDoubleArrayElements (env,B, NULL);
    cElems = (*env)-> GetDoubleArrayElements (env,C, NULL);
    assert(aElems && bElems && cElems);   

    gsl_matrix_view tempb=gsl_matrix_view_array(bElems,m,n);
    gsl_matrix_view tempc=gsl_matrix_view_array(cElems,m,n); 

    if (Uplo == jniUpper){ uplo = jniUpper;}            
    else if(Uplo == jniLower){ uplo = jniLower;}
    else{fprintf(stderr, "** Illegal Uplo setting \n"); return;}
    
    if (Side == jniLeft){
        side = jniLeft;
        lda = m;}        // C := alpha*A*B + beta*C
    else if(Side == jniRight){
        side = jniRight;
        lda = n;}    //
    else{fprintf(stderr, "** Illegal Side setting \n"); return;}
    gsl_matrix_view tempa= gsl_matrix_view_array( aElems, lda, lda); 
        
    gsl_blas_dsymm(side, uplo, alpha, &tempa.matrix, &tempb.matrix, beta, &tempc.matrix);

    (*env)-> ReleaseDoubleArrayElements (env, C, cElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, A, aElems, JNI_ABORT);
    (*env)-> ReleaseDoubleArrayElements (env, B, bElems, JNI_ABORT);   
}


JNIEXPORT void Java_JGSL_Matrix_dtrsm
(JNIEnv *env, jclass klass, jint Side, jint Uplo, jint TransA,
 jint Diag, jint m, jint n, jdouble alpha, jdoubleArray  A, jdoubleArray B){
    
    /* dtrsm: performs one of the matrix-matrix operations
     B := alpha*op( A )*B,   or   B := alpha*B*op( A )
     where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
     non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
     op( A ) = A   or   op( A ) = A**T.*/
    
    double *aElems, *bElems;
    
    int lda, i, j;
    int side, uplo, TA, diag;
    aElems = (*env)-> GetDoubleArrayElements (env,A, NULL);
    bElems = (*env)-> GetDoubleArrayElements (env,B, NULL);
    assert(aElems && bElems);    
    gsl_matrix_view tempb=gsl_matrix_view_array(bElems,m,n);

    if (TransA == jniNoTrans) {TA = jniNoTrans;}
    else if (TransA == jniTrans) {TA = jniTrans;}
    else if (TransA == jniConjTrans) {TA = jniConjTrans;}
    else {fprintf(stderr, "** Illegal TransA setting \n"); return;}
    
    if (Diag == jniNonUnit) {diag = jniNonUnit;}
    else if (Diag == jniUnit) {diag = jniUnit;}
    else {fprintf(stderr, "** Illegal Diag setting \n"); return;}

    if (Uplo == jniUpper){ uplo = jniUpper;}                //
    else if(Uplo == jniLower){ uplo = jniLower;}
    else{fprintf(stderr, "** Illegal Uplo setting \n"); return;}    

    if (Side == jniLeft){
        side = jniLeft;
        lda = m;}        // C := alpha*A*B + beta*C
    else if(Side == jniRight){
        side = jniRight;
        lda = n;}    //
    else{fprintf(stderr, "** Illegal Side setting \n"); return;}

    gsl_matrix_view tempa= gsl_matrix_view_array(aElems,lda,lda);

        
    gsl_blas_dtrsm(side, uplo, TA, diag, alpha, &tempa.matrix, &tempb.matrix);
        

    (*env)-> ReleaseDoubleArrayElements (env, B, bElems, 0);
    (*env)-> ReleaseDoubleArrayElements (env, A, aElems, JNI_ABORT);
    
}

 JNIEXPORT void JNICALL Java_JGSL_Matrix_transposeMemcpy (JNIEnv *env, jclass obj, jint m, jint n, jdoubleArray jy, jdoubleArray jx){
//	int k = 1;
	jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
 	jdouble *y = (*env)->GetDoubleArrayElements(env,jy,NULL);
 	assert(x && y);
 
 	gsl_matrix_view tempx=gsl_matrix_view_array(x,m,n);
 	gsl_matrix_view tempy=gsl_matrix_view_array(y,n,m);
 
	gsl_matrix_transpose_memcpy(&tempy.matrix, &tempx.matrix);
  
 	(*env)->ReleaseDoubleArrayElements(env,jx, x,JNI_ABORT);
 	(*env)->ReleaseDoubleArrayElements(env,jy, y, 0);

}	 

JNIEXPORT void JNICALL Java_JGSL_Matrix_transpose (JNIEnv *env, jclass obj, jint m, jdoubleArray jx){
        jdouble *x = (*env)->GetDoubleArrayElements(env,jx,NULL);
        assert(x);

        gsl_matrix_view tempx=gsl_matrix_view_array(x,m,m);
        gsl_matrix_transpose(&tempx.matrix);

        (*env)->ReleaseDoubleArrayElements(env, jx, x, 0);

}

// JNIEXPORT int JNICALL Java_JGSL_

