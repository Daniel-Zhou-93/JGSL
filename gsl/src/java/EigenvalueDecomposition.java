package JGSL;

import java.lang.*;

/** Eigenvalues and eigenvectors of a real matrix.
 <P>
    If A is symmetric, then A = V*D*V' where the eigenvalue matrix D is
    diagonal and the eigenvector matrix V is orthogonal.
    I.e. A = V.times(D.times(V.transpose())) and
    V.times(V.transpose()) equals the identity matrix.
 <P>
    If A is not symmetric, then the eigenvalue matrix D is block diagonal
    with the real eigenvalues in 1-by-1 blocks and any complex eigenvalues,
    lambda + i*mu, in 2-by-2 blocks, [lambda, mu; -mu, lambda].  The
    columns of V represent the eigenvectors in the sense that A*V = V*D,
    i.e. A.times(V) equals V.times(D).  The matrix V may be badly
    conditioned, or even singular, so the validity of the equation
    A = V*D*inverse(V) depends upon V.cond().
 **/

public class EigenvalueDecomposition  implements java.io.Serializable {
 static {
    /* load library (which will contain blas functions.)*/
    System.loadLibrary("gsl_eig");
     
 }

   /* ------------------------
    * Class variables
    * ------------------------ */
    /** Row and column dimension (square matrix).*/
    private int n;

    /** Symmetry flag.*/
    private boolean issymmetric;

    /** Arrays for internal storage of eigenvalues.
     @serial internal storage of eigenvalues.
     */
    private double[] wr, wi;
    
    /**Array for elements storage of A */
    private double[] a;
    
    /** Array for internal storage of eigenvectors.
     @serial internal storage of eigenvectors.
     */
    private double[] VR; //vector real
    private double[] VI; //vector image
    

    /* ------------------------
     * Constructor
     * ------------------------ */
    
    /** Check for symmetry, then construct the eigenvalue decomposition
     Structure to access D and V.
     @param Arg    Square matrix
     */
   
    public EigenvalueDecomposition (Matrix Arg) {
        double[][] A = Arg.getArray();
        a = Arg.getRowPackedCopy();
        n = Arg.getRowDimension();
        VR = new double[n * n];
        VI = new double[n * n];
        wr = new double[n];
        wi = new double[n];
        
        issymmetric = true;
        for (int j = 0; (j < n) & issymmetric; j++) {
            for (int i = 0; (i < n) & issymmetric; i++) {
                issymmetric = (A[i][j] == A[j][i]);
            }
        }
        
        int matrix_layout =  EigenvalueDecomposition.LAYOUT.RowMajor;
        if (issymmetric) {
            /** Calculation */
            symmv( n, a, wr, VR);
        } else{
            /** Calculation */
            nonsymmv( n, a, wr, wi, VR, VI);
        }
    }

    /* ------------------------
     * Public Methods
     * ------------------------ */
    /** Return the real eigenvector matrix */
    public Matrix getV () {
        if (issymmetric){
            return new Matrix(VR, n, n);
        } else{
            return new Matrix(VR, n, n);
        }
        
    }
    public Matrix getRealVectors () {
        if (issymmetric){
            return new Matrix(VR, n, n);
        } else{
            return new Matrix(VR, n, n);
        }
    }
    public Matrix getImagVectors(){
        if (issymmetric){
            throw new IllegalArgumentException("Matrix must nonsymmetrix.");
        }else{
            return new Matrix(VI, n, n);
        }
    }
    
    /** Return the real parts of the eigenvalues
     @return     real(diag(D))
     */
    
    public double[] getRealEigenvalues () {
        return wr;
    }
    
    /** Return the imaginary parts of the eigenvalues
     @return     imag(diag(D))
     */
    
    public double[] getImagEigenvalues () {
        return wi;
    }
    
    /** Return the block diagonal eigenvalue matrix
     @return     D
     */
    public Matrix getD () {
        Matrix X = new Matrix(n,n);
        double[][] D = X.getArray();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                D[i][j] = 0.0;
            }
            D[i][i] = wr[i];
            if (wi[i] > 0) {
                D[i][i+1] = wi[i];
            } else if (wi[i] < 0) {
                D[i][i-1] = wi[i];
            }
        }
        return X;
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
        public final static char Lower= 'L';            /** Lower triangular matrix*/
    }
    
    
    
    
    /* Eigenvector and SVD */
    public static native void nonsymmv(int n, double[] a, double[] real, double[] imag,
                                       double[] Mreal, double[] Mimag);
    
    /**inform java virtual machine that function is defined externally*/
    
    public static native void symmv( int n, double[] a, double[] eval, double[] evec);
 
    private static final long serialVersionUID = 1;
    
}







