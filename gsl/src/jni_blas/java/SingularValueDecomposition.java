package JGSL;

import java.lang.*;

/** Singular Value Decomposition.
 <P>
 For an m-by-n matrix A, the singular value decomposition is
 an m-by-m orthogonal matrix U, an m-by-n matrix which is zero except for its
 min(m,n) diagonal elements, S, and an n-by-n orthogonal matrix V so that A = U*S*V'.
 <P>
 The singular values, sigma[k] = S[k][k], are ordered so that
 sigma[0] >= sigma[1] >= ... >= sigma[n-1].
 <P>
 The singular value decompostion always exists, so the constructor will
 never fail.  The matrix condition number and the effective numerical
 rank can be computed from this decomposition.
 */

public class SingularValueDecomposition implements java.io.Serializable {
 static {
    /* load library (which will contain blas functions.)*/
    System.loadLibrary("gsl_svd");
 }

 /* ------------------------
  * Class variables
  * ------------------------ */

    /** Arrays for internal storage of U and V.*/
    private double[] a,u,v,A;

    /** Array for internal storage of singular values.*/
    private double[] s;

    /** Row and column dimensions, and min(m, n).*/
    private int m, n, l;
    
 /* ------------------------
  * Constructor
  * ------------------------ */
    public SingularValueDecomposition (Matrix Arg) {

        m = Arg.getRowDimension();

        n = Arg.getColumnDimension();
        a = new double[m * n];
        a = Arg.getRowPackedCopy();
        l = Math.min(m,n);
        s = new double[l]; //singular value
        u = new double[n * n];
        v = new double[n * n];
        int matrix_layout = SingularValueDecomposition.LAYOUT.RowMajor;
        double[] superb = new double[l];

        /** Calculation */
        dcomp(m, n, l, a, v, s);
    }

    /* ------------------------
     * Public Methods
     * ------------------------ */
    /** Return the left singular vectors U*/
    public Matrix getU () {
        Matrix X = new Matrix(m,n);
        double[][] UU = X.getArray();
        for (int j = 0; j < m; j++) {
            for (int i = 0; i < n; i++) {
                UU[j][i] = this.a[i+j*n];
            }
        }
        return X;
    }


    /** Return the right singular vectors V (V)*/
     public Matrix getV () {
         Matrix X = new Matrix(n,n);
         double[][] V = X.getArray();
         for (int j = 0; j < n; j++) {
             for (int i = 0; i < n; i++) {
                 V[j][i] = this.v[i + j*n];
             }
         }
         return X;
     }

    /** Return the one-dimensional array of singular values
     @return     diagonal of S.
     */
    
    public double[] getSingularValues () {
        return s;
    }

    /** Return the n-by-n matrix of singular values
     @return     S
     */
    public Matrix getS () {
        Matrix X = new Matrix(n,n);
        double[][] SS = X.getArray();
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < n; i++) {
                SS[i][j] = 0.0;
                SS[j][j] = this.s[j];
                }
            
        }
        return X;
    }
    public Matrix solve (Matrix B) {
        if (B.getRowDimension() != m) {
            throw new IllegalArgumentException("Matrix row dimensions must agree.");
        }
        
        // Copy right hand side
        int nx = B.getRowDimension();
        double[] b = B.getRowPackedCopy();
        double[] x = new double[B.getRowDimension()]; 
        slve( m, n, l, a, v, s, b, x);
        
        Matrix C = new Matrix(x, B.getRowDimension());
        return C;
    }

    /** Two norm
     @return     max(S)
     */
    
    public double norm2 () {
        return s[0];
    }
    
    /** Two norm condition number
     @return     max(S)/min(S)
     */
    
    public double cond () {
        return s[0]/s[l-1];
    }
    
    public int rank () {
        double eps = Math.pow(2.0,-52.0);
        double tol = Math.max(m,n)*s[0]*eps;
        int r = 0;
        for (int i = 0; i < s.length; i++) {
            if (s[i] > tol) {
                r++;
            }
        }
        return r;
    }
    private static final long serialVersionUID = 1;

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
        public final static char Lower= 'L';            /** Lower triangular matrix*/
    }
    
    
    /* Eigenvector and SVD */
    
    public static native void dcomp(int m, int n, int l, double[] a, double[] v,
                                    double[] s);
    
    public static native void slve(int m, int n, int l, double[] a, double[] v, double[] s,
                                    double[] b, double[] x);
    
    /**inform java virtual machine that function is defined externally*/
    
 }




