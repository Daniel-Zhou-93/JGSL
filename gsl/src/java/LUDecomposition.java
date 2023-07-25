package JGSL;

import java.lang.*;

public class LUDecomposition implements java.io.Serializable {
 static {
    /* load library (which will contain blas functions.)*/
//    System.loadLibrary("gsl_LUDecomposition");
    System.loadLibrary("gsl_lu");
 }
   
    /* ------------------------
     * Class variables
     * ------------------------ */
    /** Array for internal storage of decomposition.*/
    private double[] LU;

    /** Row and column dimensions, min(m, n) and INFO.*/
    private int n, m,  signum;
    private long[] p;
    
    //private double pivsign;

    /** Internal storage of pivot vector.*/
    //private int[] ipiv, piv;

    /* ------------------------
     * Constructor
     * ------------------------ */
    public LUDecomposition (Matrix Arg) {
        
        LU = Arg.getRowPackedCopy();
        m = Arg.getRowDimension();
        n = Arg.getColumnDimension();
        //l = Math.min(m, n);
        signum = 1;
        p = new long[n];
        //int temp;
        
        //for (int i = 0; i < l; i++) {
        //    ipiv[i] = i;
        //    piv[i] = i;
        //}
        //int matrix_layout = LUDecomposition.LAYOUT.RowMajor;
        signum= decomp( n, LU, p,  signum);

        /*for (int i = 0; i < l; i++){
            temp = piv[ipiv[i] - 1];
            piv[ipiv[i] - 1] = piv[i];
            piv[i] = temp;
            if (piv[i] != i ){pivsign = -pivsign;}
        }
	}*/

    /* ------------------------
     * Public Methods
     * ------------------------ */

    /** Is the matrix nonsingular?*/
    /*public boolean isNonsingular () {
        if(info == 0 ){
            return true;
        } else if (info > 0){
            return false;
        } else {
            return false;
        }
        /* should I add check m == n? */
    }

   /** Return lower triangular factor */
    public Matrix getL () {
        //int nrow = m, ncol = ((m >= n) ? n: m);
        Matrix X = new Matrix(n, n);
        double[][] L = X.getArray();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (j > i) {
                    L[j][i] = LU[i + j * n];
                } else if (i == j) {
                    L[j][i] = 1.0;
                } else {
                    L[j][i] = 0.0;
                }
            }
        }
        return X;
    }

    /** Return upper triangular factor*/
    
    public Matrix getU () {
        //int nrow = Math.min(m, n), ncol = n;
        Matrix X = new Matrix(n, n);
        double[][] U = X.getArray();
        
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++) {
                U[j][i] = ((j <= i) ? LU[i + j * n] : 0.0);
            }
        }
        return X;
    }
        
    /** Return pivot permutation vector */
    public long[] getPivot () {
        return p;
    }
    
    /** Return pivot permutation vector as a one-dimensional double array
     @return     (double) piv
     */
    
    /*public double[] getDoublePivot () {
        double[] vals = new double[l];
        for (int i = 0; i < l; i++) {
            vals[i] = (double) piv[i];
        }
        return vals;
    }*/

   /** Return inverse of matrix */
    public Matrix inver () {
        //int nrow = m, ncol = ((m >= n) ? n: m);

        double[] b = new double[n*n] ;// B is a vector
        invert(n, LU, b, p);
        Matrix X = new Matrix(b, n, n);
        return X;
    }


    /** Determinant
     @return     det(A)
     @exception  IllegalArgumentException  Matrix must be square
     */

    public double det( ) {
        if (m != n) {
            throw new IllegalArgumentException("Matrix must be square.");
        }
        
	return(det(n, LU, signum));
    }
    
    public double lndet( ) {
        if (m != n) {
            throw new IllegalArgumentException("Matrix must be square.");
        }
        
	return(lndet(n, LU));
    }
    
    /** Solve A * X = B */
    public Matrix solve (Matrix B) {
        //if (B.getRowDimension() != m) {
        //    throw new IllegalArgumentException("Matrix row dimensions must agree.");
        //}
        //if (m != n) {
        //    throw new IllegalArgumentException("Matrix must be square.");
        //}
        /*if (!this.isNonsingular()) {
            throw new RuntimeException("Matrix is singular.");
        }*/
        
        //int matrix_layout = LUDecomposition.LAYOUT.ColMajor;
        //int Trans = LUDecomposition.TRANSPOSE.NoTrans;
        //int nrhs = B.getColumnDimension();
        //int ldb = m;
        //int lda = m;
        double[] b = B.getRowPackedCopy();// B is a vector
        double[] x = new double[n];
        slve(n, LU, b, p, x);
        Matrix C = new Matrix(x, n);
        return C;
    }

  

    public final static class LAYOUT {
        private LAYOUT() {}
        public final static int RowMajor = 101;
        public final static int ColMajor = 102;
    }
    
    public final static class TRANSPOSE {
        private TRANSPOSE() {}
        public final static int NoTrans = 111;         /** trans='N' */
        public final static int Trans= 112;            /** trans='T' */
        public final static int ConjTrans= 113;        /** trans='C' */
    }
    
    public final static class UPLO {
        private UPLO() {}
        public final static int Upper = 121;           /** Upper triangular matrix */
        public final static int Lower= 122;            /** Lower triangular matrix*/
    }
    
    
    /* LU */
    public static native int decomp( int n, double[] a, long[] p, int signum);
    
    public static native int slve( int n, double[] a, double[] b, long[] p, double[] x);

    public static native double det(int n, double[] a, int signum);

    public static native double lndet(int n, double[] a);

    public static native int invert( int n, double[] a, double[] b, long[] p);
    
    private static final long serialVersionUID = 1;
}

