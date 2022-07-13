package JGSL;

import java.text.NumberFormat;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;
import java.text.FieldPosition;
import java.io.PrintWriter;
import java.io.BufferedReader;
import java.io.StreamTokenizer;


/**BLAS.java*/

public class Matrix implements Cloneable, java.io.Serializable {
 private Matrix() {}
    static {
    /* load library (which will contain blas functions.)*/
    System.loadLibrary("gsl_blas");
 }
    
    /**inform java virtual machine that function is defined externally*/
    
    /* -----------------------
     * Class variables
     * ---------------------- */
    
    private double[][] A;
    
    private int m, n;
    
    /* -----------------------
     * Constructors
     * ----------------------- */
    
    /** Construct an m-by-n matrix of zeros. */
    public Matrix (int m, int n) {
        this.m = m;
        this.n = n;
        A = new double[m][n];
    }
    
    /** Construct an m-by-n constant matrix. */
    public Matrix (int m, int n, double s) {
        this.m = m;
        this.n = n;
        A = new double[m][n];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                A[i][j] = s;
            }
        }
    }
    
    /** Construct a matrix from a 2-D array. */
    public Matrix (double[][] A) {
        m = A.length;
        n = A[0].length;
        for (int i = 0; i < m; i++) {
            if (A[i].length != n) {
                throw new IllegalArgumentException("All rows must have the same length.");
            }
        }
        this.A = A;
    }
    
    /** Construct a matrix quickly without checking arguments. */
    public Matrix (double[][] A, int m, int n) {
        this.A = A;
        this.m = m;
        this.n = n;
    }
    
    /** Construct a matrix from a one-dimensional packed array
     @param vals One-dimensional array of doubles, packed by columns (ala Fortran).
     @param m    Number of rows.
     @exception  IllegalArgumentException Array length must be a multiple of m.
     */
    
    public Matrix (double vals[], int m) {
        this.m = m;
        n = (m != 0 ? vals.length/m : 0);
        if (m*n != vals.length) {
            throw new IllegalArgumentException("Array length must be a multiple of m.");
        }
        A = new double[m][n];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                A[i][j] = vals[i+j*m];
            }
        }
    }

    public Matrix (double vals[], int m, int n) {
        this.m = m;
        this.n = n;
        A = new double[m][n];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                A[i][j] = vals[i*n+j];
            }
        }
    }
    
    /* ------------------------
     *    Public Methods
     * ------------------------ */
    
    /** Construct a matrix from a copy of a 2-D array.*/
    
    public static Matrix constructWithCopy(double[][] A) {
        int m = A.length;
        int n = A[0].length;
        Matrix X = new Matrix(m,n);
        double[][] C = X.getArray();
        for (int i = 0; i < m; i++) {
            if (A[i].length != n) {
                throw new IllegalArgumentException
                ("All rows must have the same length.");
            }
            for (int j = 0; j < n; j++) {
                C[i][j] = A[i][j];
            }
        }
        return X;
    }
    
    
    /** Make a deep copy of a matrix
     */
    
    public Matrix copy () {
        Matrix X = new Matrix(m,n);
        double[][] C = X.getArray();
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                C[i][j] = A[i][j];
            }
        }
        return X;
    }
    
    /** Clone the Matrix object.
     */
    
    public Object clone () {
        return this.copy();
    }
    
    /** Access the internal two-dimensional array. */
    public double[][] getArray () {
        return A;
    }

     /** Copy the internal two-dimensional array.*/
     public double[][] getArrayCopy () {
	 double[][] C = new double[m][n];
         for (int i = 0; i < m; i++) {	 
		 for (int j = 0; j < n; j++) {
			 C[i][j] = A[i][j];
		 }
	 }
	 return C;
     }

    /** Make a one-dimensional column packed copy of the internal array.*/
    public double[] getColumnPackedCopy () {
        double[] vals = new double[m*n];
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < m; i++) {
                vals[i + j*m] = A[i][j];
            }
        }
        return vals;
    }
    
    /** Make a one-dimensional row packed copy of the internal array. */
    public double[] getRowPackedCopy() {
        double[] vals = new double[m*n];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                vals[i*n + j] = A[i][j];
            }
        }
        return vals;
    }

    /** Get row dimension.*/
    public int getRowDimension () {
        return m;
    }
    /** Get column dimension.*/
    public int getColumnDimension () {
        return n;
    }
    /** Get a single element.*/
    public double get (int i, int j) {
        return A[i][j];
    }
    
    public Matrix getMatrix (int i0, int i1, int j0, int j1) {
        Matrix X = new Matrix(i1-i0+1,j1-j0+1);
        double[][] B = X.getArray();
        try {
            for (int i = i0; i <= i1; i++) {
                for (int j = j0; j <= j1; j++) {
                    B[i-i0][j-j0] = A[i][j];
                }
            }
        } catch(ArrayIndexOutOfBoundsException e) {
            throw new ArrayIndexOutOfBoundsException("Submatrix indices");
        }
        return X;
    }

    
    /** Get a submatrix.
     @param r    Array of row indices.
     @param c    Array of column indices.
     @return     A(r(:),c(:))
     @exception  ArrayIndexOutOfBoundsException Submatrix indices
     */
    
    public Matrix getMatrix (int[] r, int[] c) {
        Matrix X = new Matrix(r.length,c.length);
        double[][] B = X.getArray();
        try {
            for (int i = 0; i < r.length; i++) {
                for (int j = 0; j < c.length; j++) {
                    B[i][j] = A[r[i]][c[j]];
                }
            }
        } catch(ArrayIndexOutOfBoundsException e) {
            throw new ArrayIndexOutOfBoundsException("Submatrix indices");
        }
        return X;
    }


    
   
    /** Get a submatrix.
     @param i0   Initial row index
     @param i1   Final row index
     @param c    Array of column indices.
     @return     A(i0:i1,c(:))
     @exception  ArrayIndexOutOfBoundsException Submatrix indices
     */
    
    public Matrix getMatrix (int i0, int i1, int[] c) {
        Matrix X = new Matrix(i1-i0+1,c.length);
        double[][] B = X.getArray();
        try {
            for (int i = i0; i <= i1; i++) {
                for (int j = 0; j < c.length; j++) {
                    B[i-i0][j] = A[i][c[j]];
                }
            }
        } catch(ArrayIndexOutOfBoundsException e) {
            throw new ArrayIndexOutOfBoundsException("Submatrix indices");
        }
        return X;
    }

    /** Get a submatrix.
     @param r    Array of row indices.
     @param j0   Initial column index
     @param j1   Final column index
     @return     A(r(:),j0:j1)
     @exception  ArrayIndexOutOfBoundsException Submatrix indices
     */
    
    public Matrix getMatrix (int[] r, int j0, int j1) {
        Matrix X = new Matrix(r.length,j1-j0+1);
        double[][] B = X.getArray();
        try {
            for (int i = 0; i < r.length; i++) {
                for (int j = j0; j <= j1; j++) {
                    B[i][j-j0] = A[r[i]][j];
                }
            }
        } catch(ArrayIndexOutOfBoundsException e) {
            throw new ArrayIndexOutOfBoundsException("Submatrix indices");
        }
        return X;
    }

    public Matrix getMatrix (long[] r, int j0, int j1) {
        Matrix X = new Matrix(r.length,j1-j0+1);
        double[][] B = X.getArray();
        try {
            for (int i = 0; i < r.length; i++) {
                for (int j = j0; j <= j1; j++) {
                    B[i][j-j0] = A[(int)r[i]][j];
                }
            }
        } catch(ArrayIndexOutOfBoundsException e) {
            throw new ArrayIndexOutOfBoundsException("Submatrix indices");
        }
        return X;
    }
    
    /** Set a single element.
     @param i    Row index.
     @param j    Column index.
     @param s    A(i,j).
     @exception  ArrayIndexOutOfBoundsException
     */
    
    public void set (int i, int j, double s) {
        A[i][j] = s;
    }
    
    /** Set a submatrix.
     @param i0   Initial row index
     @param i1   Final row index
     @param j0   Initial column index
     @param j1   Final column index
     @param X    A(i0:i1,j0:j1)
     @exception  ArrayIndexOutOfBoundsException Submatrix indices
     */
    
    public void setMatrix (int i0, int i1, int j0, int j1, Matrix X) {
        try {
            for (int i = i0; i <= i1; i++) {
                for (int j = j0; j <= j1; j++) {
                    A[i][j] = X.get(i-i0,j-j0);
                }
            }
        } catch(ArrayIndexOutOfBoundsException e) {
            throw new ArrayIndexOutOfBoundsException("Submatrix indices");
        }
    }
    
    /** Set a submatrix.
     @param r    Array of row indices.
     @param c    Array of column indices.
     @param X    A(r(:),c(:))
     @exception  ArrayIndexOutOfBoundsException Submatrix indices
     */
    
    public void setMatrix (int[] r, int[] c, Matrix X) {
        try {
            for (int i = 0; i < r.length; i++) {
                for (int j = 0; j < c.length; j++) {
                    A[r[i]][c[j]] = X.get(i,j);
                }
            }
        } catch(ArrayIndexOutOfBoundsException e) {
            throw new ArrayIndexOutOfBoundsException("Submatrix indices");
        }
    }
    
    /** Set a submatrix.
     @param r    Array of row indices.
     @param j0   Initial column index
     @param j1   Final column index
     @param X    A(r(:),j0:j1)
     @exception  ArrayIndexOutOfBoundsException Submatrix indices
     */
    
    public void setMatrix (int[] r, int j0, int j1, Matrix X) {
        try {
            for (int i = 0; i < r.length; i++) {
                for (int j = j0; j <= j1; j++) {
                    A[r[i]][j] = X.get(i,j-j0);
                }
            }
        } catch(ArrayIndexOutOfBoundsException e) {
            throw new ArrayIndexOutOfBoundsException("Submatrix indices");
        }
    }
    
    /** Set a submatrix.
     @param i0   Initial row index
     @param i1   Final row index
     @param c    Array of column indices.
     @param X    A(i0:i1,c(:))
     @exception  ArrayIndexOutOfBoundsException Submatrix indices
     */
    
    public void setMatrix (int i0, int i1, int[] c, Matrix X) {
        try {
            for (int i = i0; i <= i1; i++) {
                for (int j = 0; j < c.length; j++) {
                    A[i][c[j]] = X.get(i-i0,j);
                }
            }
        } catch(ArrayIndexOutOfBoundsException e) {
            throw new ArrayIndexOutOfBoundsException("Submatrix indices");
        }
    }
    
    /** Matrix transpose.
     @return    A'
     
    
    public Matrix transpose () {
        Matrix X = new Matrix(n,m);
        double[][] C = X.getArray();
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                C[j][i] = A[i][j];
            }
        }
        return X;
    }
   
	**/

    public Matrix TransposeMemcpy () {    
	double[] x = this.getRowPackedCopy();
	double[] y = new double[n*m];
 	transposeMemcpy(m, n, y, x);
	Matrix Y = new Matrix(y, n, m);
	return Y;
    }

    public Matrix Transpose () {
	    if (this.m != this.n) {
		throw new ArrayIndexOutOfBoundsException("This function needs square matrices");
	    }	
        double[] x = this.getRowPackedCopy();
        transpose(m, x);
        Matrix X = new Matrix(x, n, m);
        return X;
    }

    
    /**  Unary minus
     @return    -A
     */
    
    public Matrix uminus ( ) {
        double[] a = this.getRowPackedCopy();
        dscal(m * n, -1, a);
        Matrix C = new Matrix(a, m, n);
        return C;
    }
    
    /* No function in Level 1-3 involves operation of addition, operation of matrix
     addition is written without using functions in Level 1-3        */
    /* C=A+B     */
    public Matrix plus (Matrix B){
        checkMatrixDimensions(B);
        double[] a = this.getRowPackedCopy();
        double[] b = B.getRowPackedCopy();
        daxpy(m * n, 1, b, a);
        Matrix C = new Matrix (a , m, n);
        return C;
    }
    
    
    /** A = A + B
     @param B    another matrix
     @return     A + B
     */
    
    public Matrix plusEquals (Matrix B) {
        checkMatrixDimensions(B);
        double[] a = this.getRowPackedCopy();
        double[] b = B.getRowPackedCopy();
        daxpy(m * n, 1, b, a);
	for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                A[i][j] = a[i*n+j];
            }
        }
	return this;
    }

    /* C = A - B     */
    public Matrix minus (Matrix B){
        checkMatrixDimensions(B);
        double[] a = this.getRowPackedCopy();
        double[] b = B.getRowPackedCopy();
        daxpy(m * n, -1.0, b, a);
        Matrix X = new Matrix (a , m, n);
        return X;
    }
    
    /** A = A - B
     @param B    another matrix
     @return     A - B
     */
 
    
    public Matrix minusEquals (Matrix B) {
        checkMatrixDimensions(B);
        double[] a = this.getRowPackedCopy();
        double[] b = B.getRowPackedCopy();
        daxpy(m * n, -1.0, b, a);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                A[i][j] = a[i*n+j];
            }
        }
	return this;
    }

    
    
    /*Declaration
     B = alpha * A
     */
    
    public Matrix times (double alpha){
        double[] a = this.getRowPackedCopy();
        dscal(m * n, alpha, a);
        Matrix X = new Matrix(a, m, n);
        return X;
    }
    
    
    /*Declaration
     X = A * B
     */
    
    public Matrix times (Matrix B) {
        if (B.m != n) {
            throw new IllegalArgumentException("Matrix dimensions must agree.");
        }
        double[] a = this.getRowPackedCopy();
        double[] b = B.getRowPackedCopy();
        double[] c = new double[m * B.getColumnDimension()];
        dgemm(Matrix.TRANSPOSE.NoTrans, Matrix.TRANSPOSE.NoTrans,
		m, n, B.getRowDimension(), B.getColumnDimension(), 1, a, 			b, 0, c);
        Matrix X = new Matrix(c, m, B.getColumnDimension());
        return X;
    }

    public Matrix times (Matrix B, int TransA, int TransB) {
        int nrow, ncol, bnrow, bncol;
        
        if (TransA == Matrix.TRANSPOSE.NoTrans){
            nrow = m;
            ncol = n;
        }
        else{
            nrow = n;
            ncol = m;
        }
        
        if (TransB == Matrix.TRANSPOSE.NoTrans){
            bnrow = B.m;
            bncol = B.n;
        }
        else{
            bnrow = B.n;
            bncol = B.m;
        }

        if (ncol != bnrow) {
            throw new IllegalArgumentException("Matrix dimensions must agree.");
        }

        double[] a = this.getRowPackedCopy();
        double[] b = B.getRowPackedCopy();
        double[] c = new double[nrow * bncol];
        dgemm( TransA, TransB,
              m, n, B.m, B.n, 1, a, b, 0, c);
        Matrix X = new Matrix(c, nrow, bncol);
        return X;
    }
    
    
    
    
    /* Using daxpy, constants times a matrix plus a matrix
     C = A + alpha * B
     */
    public Matrix plus (Matrix B, double alpha){
        checkMatrixDimensions(B);
        double[] a = this.getRowPackedCopy();
        double[] b = B.getRowPackedCopy();
        daxpy(m * n, alpha, b, a);
        Matrix C = new Matrix (a , m, n);
        return C;
    }
    
    
    /* Using dtrmm, C=alpha*op(A)*B or  C=alpha*B*op(A)        */
    public  Matrix tritimes (Matrix B, double alpha){
        checkMatrixDimensions(B);
        double[] a = this.getRowPackedCopy();
        double[] b = B.getRowPackedCopy();
        /* need to check if A is square matrix*/
        dtrmm( Matrix.SIDE.Left, Matrix.UPLO.Upper, Matrix.TRANSPOSE.NoTrans, Matrix.DIAG.NonUnit, B.getRowDimension(), B.getColumnDimension(), alpha, a, b);
        Matrix C = new Matrix (b, B.getRowDimension(), B.getColumnDimension());
        return C;
    }
    
    
    /** Generate identity matrix
     @param m    Number of rows.
     @param n    Number of colums.
     @return     An m-by-n matrix with ones on the diagonal and zeros elsewhere.
     */
    
    public static Matrix identity (int m, int n) {
        Matrix A = new Matrix(m,n);
        double[][] X = A.getArray();
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                X[i][j] = (i == j ? 1.0 : 0.0);
            }
        }
        return A;
    }
   

   /** Print the matrix to stdout - row-wise.
    *Really simple - unformatted output  
   */
    public void printmat() {
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++)
                System.out.print(A[i][j]+" ");
            System.out.print("\n");
        }
            System.out.println();
    }
 
    
    /** Print the matrix to stdout.   Line the elements up in columns
     * with a Fortran-like 'Fw.d' style format.
     @param w    Column width.
     @param d    Number of digits after the decimal.
     */
    
    public void print (int w, int d) {
        print(new PrintWriter(System.out,true),w,d); }
    
    /** Print the matrix to the output stream.   Line the elements up in
     * columns with a Fortran-like 'Fw.d' style format.
     @param output Output stream.
     @param w      Column width.
     @param d      Number of digits after the decimal.
     */
    
    public void print (PrintWriter output, int w, int d) {
        DecimalFormat format = new DecimalFormat();
        format.setDecimalFormatSymbols(new DecimalFormatSymbols(Locale.US));
        format.setMinimumIntegerDigits(1);
        format.setMaximumFractionDigits(d);
        format.setMinimumFractionDigits(d);
        format.setGroupingUsed(false);
        print(output,format,w+2);
    }
    
    /** Print the matrix to stdout.  Line the elements up in columns.
     * Use the format object, and right justify within columns of width
     * characters.
     * Note that is the matrix is to be read back in, you probably will want
     * to use a NumberFormat that is set to US Locale.
     @param format A  Formatting object for individual elements.
     @param width     Field width for each column.
     @see java.text.DecimalFormat#setDecimalFormatSymbols
     */
    
    public void print (NumberFormat format, int width) {
        print(new PrintWriter(System.out,true),format,width); }
    
    // DecimalFormat is a little disappointing coming from Fortran or C's printf.
    // Since it doesn't pad on the left, the elements will come out different
    // widths.  Consequently, we'll pass the desired column width in as an
    // argument and do the extra padding ourselves.
    
    /** Print the matrix to the output stream.  Line the elements up in columns.
     * Use the format object, and right justify within columns of width
     * characters.
     * Note that is the matrix is to be read back in, you probably will want
     * to use a NumberFormat that is set to US Locale.
     @param output the output stream.
     @param format A formatting object to format the matrix elements
     @param width  Column width.
     @see java.text.DecimalFormat#setDecimalFormatSymbols
     */
    
    public void print (PrintWriter output, NumberFormat format, int width) {
        output.println();  // start on new line.
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                String s = format.format(A[i][j]); // format the number
                int padding = Math.max(1,width-s.length()); // At _least_ 1 space
                for (int k = 0; k < padding; k++)
                    output.print(' ');
                output.print(s);
            }
            output.println();
        }
        output.println();   // end with blank line.
    }
    
    /** Read a matrix from a stream.  The format is the same the print method,
     * so printed matrices can be read back in (provided they were printed using
     * US Locale).  Elements are separated by
     * whitespace, all the elements for each row appear on a single line,
     * the last row is followed by a blank line.
     @param input the input stream.
     */
    
    public static Matrix read (BufferedReader input) throws java.io.IOException {
        StreamTokenizer tokenizer= new StreamTokenizer(input);
        
        // Although StreamTokenizer will parse numbers, it doesn't recognize
        // scientific notation (E or D); however, Double.valueOf does.
        // The strategy here is to disable StreamTokenizer's number parsing.
        // We'll only get whitespace delimited words, EOL's and EOF's.
        // These words should all be numbers, for Double.valueOf to parse.
        
        tokenizer.resetSyntax();
        tokenizer.wordChars(0,255);
        tokenizer.whitespaceChars(0, ' ');
        tokenizer.eolIsSignificant(true);
        java.util.Vector<Double> vD = new java.util.Vector<Double>();
        
        // Ignore initial empty lines
        while (tokenizer.nextToken() == StreamTokenizer.TT_EOL);
        if (tokenizer.ttype == StreamTokenizer.TT_EOF)
            throw new java.io.IOException("Unexpected EOF on matrix read.");
        do {
            vD.addElement(Double.valueOf(tokenizer.sval)); // Read & store 1st row.
        } while (tokenizer.nextToken() == StreamTokenizer.TT_WORD);
        
        int n = vD.size();  // Now we've got the number of columns!
        double row[] = new double[n];
        for (int j=0; j<n; j++)  // extract the elements of the 1st row.
            row[j]=vD.elementAt(j).doubleValue();
        java.util.Vector<double[]> v = new java.util.Vector<double[]>();
        v.addElement(row);  // Start storing rows instead of columns.
        while (tokenizer.nextToken() == StreamTokenizer.TT_WORD) {
            // While non-empty lines
            v.addElement(row = new double[n]);
            int j = 0;
            do {
                if (j >= n) throw new java.io.IOException
                    ("Row " + v.size() + " is too long.");
                row[j++] = Double.valueOf(tokenizer.sval).doubleValue();
            } while (tokenizer.nextToken() == StreamTokenizer.TT_WORD);
            if (j < n) throw new java.io.IOException
                ("Row " + v.size() + " is too short.");
        }
        int m = v.size();  // Now we've got the number of rows.
        double[][] A = new double[m][];
        v.copyInto(A);  // copy the rows out of the vector
        return new Matrix(A);
    }

    
    /* ------------------------
     Private Methods
     * ------------------------ */
    
    /** Check if size(A) == size(B) **/
    
    private void checkMatrixDimensions (Matrix B) {
        if (B.m != m || B.n != n) {
            throw new IllegalArgumentException("Matrix dimensions must agree.");
        }
    }
    
    private static final long serialVersionUID = 1;
    
    
    public static class LAYOUT {
        private LAYOUT() {}
        /** row-major arrays */
        public final static int RowMajor= 101;
        /** column-major arrays */
        public final static int ColMajor= 102;
    }
    
    public static class TRANSPOSE {
        private TRANSPOSE() {}
        /** trans = 'N' */
        public final static int NoTrans = 111;
        /** trans = 'T' */
        public final static int Trans= 112;
        /** trans = 'C'*/
        public final static int ConjTrans= 113;
    }

    public static class UPLO {
        private UPLO() {}
        /** Upper triangular matrix */
        public final static int Upper = 121;
        /** Lower triangular matrix*/
        public final static int Lower = 122;
    }
    
    public static class DIAG {
        private DIAG() {}
        /** not assumed to be unit  */
        public final static int NonUnit = 131;
        /** assumed to be unit */
        public final static int Unit = 132;
    }
    
    public static class SIDE {
        private SIDE() {}
        /** B := alpha*op( A )*B */
        public final static int Left = 141;
        /** B := alpha*B*op( A ) */
        public final static int Right = 142;
    }
    
    /* Level 1: */
    public static native void dscal( int n, double alpha, double[] x);
    
    public static native void daxpy( int n, double alpha, double[] x, double[] y);
    
    public static native double ddot(int n, double[] x, double[] y);
    
    /* Level 2: */
    public static native void dgemv(int Trans, int m, int n, double
                                    alpha, double[] A, double[] x, double beta, 
                                    double[] y);
    
    public static native void dtrmv(int Uplo, int Trans, int Diag,
                                    int n, double[] A, double[] x);
    
    public static native void dsymv(int Uplo, int n, double alpha,
                                    double[] A, double[] x, double beta, double[] y);
    
    /* Level 3: */

    public static native void dgemm(int TransA, int TransB, int m,
                                    int n, int k, int l, double alpha, double[] A,
                                    double[] B, double beta, double[] C);
    
    public static native void dtrmm(int Side, int Uplo, int TransA,
                                    int Diag, int m, int n, double alpha,
                                    double[] A, double[] B);
    
    public static native void dtrsm(int Side, int Uplo, int TransA,
                                    int Diag, int m, int n, double alpha,
                                    double[] A, double[] B);
    
    public static native void dsymm(int Side, int Uplo, int m, int n,
                                    double alpha, double[] A, double[] B,
                                    double beta, double[] C);

   public static native void transposeMemcpy(int m, int n, double[] y, double[] x);	

   public static native void transpose(int m, double[] x);	


}

