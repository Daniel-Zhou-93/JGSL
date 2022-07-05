package JGSL;

import JGSL.*;
import java.lang.*;
/** QR Decomposition.
 <P>
    For an m-by-n matrix A with m >= n, the QR decomposition is an m-by-n
    orthogonal dormqrmatrix Q and an n-by-n upper triangular matrix R so that
    A = Q*R.
 <P>
    The QR decompostion always exists, even if the matrix does not have
    full rank, so the constructor will never fail.  The primary use of the
    QR decomposition is in the least squares solution of nonsquare systems
    of simultaneous linear equations.  This will fail if isFullRank()
    returns false.
 */

public class QRDecomposition implements java.io.Serializable {
 static {
    /* load library (which will contain lapack functions.)*/
    System.loadLibrary("gsl_qr");
 }
     
    /* ------------------------
     *  Class variables
     *  ------------------------ */
    /** Array for internal storage of decomposition.*/
    private double[] QR;
    private double[] tau;
    
    /** Row and column dimensions, and min(m, n).*/
    private int m, n, l;

    /* ------------------------
     Constructor
     * ------------------------ */

    /** Check for symmetry, then construct the eigenvalue decomposition
     Structure to access D and V.
     @exception  IllegalArgumentException  Matrix row dimensions must agree.
     */
    public QRDecomposition (Matrix Arg) {
        m = Arg.getRowDimension();
        n = Arg.getColumnDimension();
        l = Math.min(m, n);

        QR = Arg.getRowPackedCopy();
        tau = new double[l];
	int matrix_layout = QRDecomposition.LAYOUT.RowMajor;
        /** Calculation */
        dcomp(m, n, l, QR, tau);
        
    }

    /* ------------------------
     * Public Methods
     * ------------------------ */

    /** Is the matrix full column rank?
     @return     true if R, and hence A, has full rank.
     */
    
    public boolean isFullRank () {
        for(int j = 0; j < l; j++){
            if (QR[j * m + j] == 0)
                return false;
        }
        return true;
    }
    
    /** no get H */
    
    /** Return the upper triangular factor*/
    public Matrix getR () {
        Matrix X = new Matrix(m,n);
	    double[][] R = X.getArray();
	    for (int j = 0; j < m; j++) {
		    for (int i = 0; i < n; i++) {
                if(j <= i){
                    R[j][i] = QR[i+j*n];
                } else {
                    R[j][i] = 0.0;
                }
		    }
	    }
	    return X;

    }


   /** Generate and return the (economy-sized) orthogonal factor*/
    public Matrix getQ () {
        int lda = m;
        int matrix_layout = QRDecomposition.LAYOUT.RowMajor;
        double[] QRtemp = new double[m * n];
        double[] Q = new double[m * m];
        for(int i = 0; i < (m * n); i++){
            QRtemp[i] = QR[i];
        }
        /** Calculation */

        unpack( m, n, l, QRtemp, tau, Q);
        return new Matrix(Q, m, m);
        /** here n has to be less than m; otherwise dorgqr will fail ?*/
    }

    /** Least squares solution of A*X = B
     @param B    A Matrix with as many rows as A and any number of columns.
     @return     X that minimizes the two norm of Q*R*X-B.
     @exception  IllegalArgumentException  Matrix row dimensions must agree.
     @exception  RuntimeException  Matrix is rank deficient.
     */

    public Matrix solve (Matrix B) {
        //if (B.getRowDimension() != m) {
        //    throw new IllegalArgumentException("Matrix row dimensions must agree.");
        //}
        //if (!this.isFullRank()) {
        //    throw new RuntimeException("Matrix is rank deficient.");
        //}
        
        // Copy right hand side
        int nx = B.getRowDimension();
        double[] b = B.getRowPackedCopy();
        double[] x = new double[B.getRowDimension()]; 
        slve( m, n, l, QR, tau, b, x);
        
        Matrix C = new Matrix(x, B.getRowDimension());
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
    
    public final static class SIDE {
        private SIDE() {}
        public final static char Left = 'L';            /** apply Q or Q**T from the Left */
        public final static char Right= 'R';            /** apply Q or Q**T from the Right */
    }
    
    public final static class UPLO {
        private UPLO() {}
        public final static char Upper = 'U';           /** Upper triangular matrix */
        public final static char Lower= 'L';            /** Lower triangular matrix*/
    }
    
    
    
    /* QR */
    public static native int dcomp( int m, int n, int l, double[] a, double[] tau);
    
    public static native void unpack( int m, int n, int l, double[] a, double[] tau, double[] q);
    
    public static native void slve( int m, int n, int l, double[] a, double[] tau,
                                    double[] b, double[] x);


    private static final long serialVersionUID = 1;
    /**inform java virtual machine that function is defined externally*/
    
}







