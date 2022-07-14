package JGSL;

import JGSL.*;
import java.lang.*;

public class CholeskyDecomposition implements java.io.Serializable {
	static {
        /* load library (which will contain blas functions.)*/
         System.loadLibrary("gsl_chol");
	}

    /* ------------------------
     *  Class variables
     *  ------------------------ */
    
    /** Array for internal storage of decomposition.
     @serial internal array storage.
     */
    private double[] a;

    /** Row and column dimension (square matrix).*/
	private int n, info=0;

    /** Symmetric and positive definite flag.*/
	private boolean isspd;

    /* ------------------------
     * Constructor
     * ------------------------ */
    
    /** Cholesky algorithm for symmetric and positive definite matrix.
     Structure to access L and isspd flag.
     @param  Arg   Square, symmetric matrix.
     */
    
	public CholeskyDecomposition (Matrix Arg) {
        
        double[][] A = Arg.getArray();
        a = Arg.getRowPackedCopy();
        n = Arg.getRowDimension();
        isspd = (Arg.getColumnDimension() == n);
        
        for (int j = 0; j < n; j++){
            for (int k = 0; k < j; k++){
                isspd = isspd & (A[k][j] == A[j][k]);
            }
        }
        
        
        int matrix_layout = CholeskyDecomposition.LAYOUT.RowMajor;
        dcomp(n, a, info);
        isspd = isspd & (info == 1);
    }

   /* ------------------------
    * Public Methods
    * ------------------------ */

    /** Is the matrix symmetric and positive definite?*/
	public boolean isSPD () {
		return isspd;
	}

	/** Return triangular factor.*/
	public Matrix getL () {
		Matrix X = new Matrix(n,n);
		double[][] L = X.getArray();
        for (int j = 0; j < n; j++) 
		for (int i = 0; i <= j; i++) {
			L[j][i] = this.a[j*n+i];   /* check whether the switch of j and i can be avoid of not*/
		}
	return X;
	}
        
	/** Solve A*X = B
     @param  B   A Matrix with as many rows as A and any number of columns.
     @return     X so that L*L'*X = B
     @exception  IllegalArgumentException  Matrix row dimensions must agree.
     @exception  RuntimeException  Matrix is not symmetric positive definite.
     */
    
	public Matrix solve (Matrix B) {
	//	if (B.getRowDimension() != n) {
	//		throw new IllegalArgumentException("Matrix row dimensions must agree.");
        //}
	//	if (!isspd) {
	//		throw new RuntimeException("Matrix is not symmetric positive definite.");
	//	}
        /*here the isspd only check square matrix...*/
        
        int matrix_layout = CholeskyDecomposition.LAYOUT.RowMajor;
	double[] b = B.getRowPackedCopy();
        double[] x = new double[n];
	slve( n, a, b, x);
	Matrix C = new Matrix(x,n);
	return C;
	}
 

    public final static class LAYOUT {
        private LAYOUT() {}
        public final static int RowMajor = 101;
        public final static int ColMajor = 102;
    }

    public final static class TRANSPOSE {
        private TRANSPOSE() {}
        public final static char NoTrans = 'N';         /** trans='N' */
        public final static char Trans= 'T';            /** trans='T' */
        public final static char ConjTrans= 'C';        /** trans='C' */
    }

    public final static class UPLO {
        private UPLO() {}
        public final static char Upper = 'U';           /** Upper triangular matrix */
        public final static char Lower = 'L';            /** Lower triangular matrix*/
    }


    private static final long serialVersionUID = 1;
    
    /**inform java virtual machine that function is defined externally*/
    public static native int dcomp(int n, double[] a, int info);
    
    public static native void invert(int n, double[] a);
    
    public static native void slve(int n, double[] a, double[] b, double[] x);
    
}




