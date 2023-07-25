import JGSL.*;
import java.io.*;
import java.util.zip.GZIPInputStream;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;


public final class JGSLExamplesLU {
	private JGSLExamplesLU() {}
	public static void main(String[] args) {
        //
        // Prepare the matrices and other parameters
        //
		System.out.println("###   Testfile for jgsl   ### \n \n");
		System.out.println("###   Parameter Preparation   ###");
		int M=3, N=3, K=3;
		double[] A = new double[] {12, -51, 4, 6, 167, -68, -4, 24, -41};
		double[] B = new double[] {4, -2, -6, -2, 10, 9, -6, 9, 14};
		double[] D = new double[9];
		double[] D2 = new double[9];
		int[] ipiv = new int[] {0, 0, 0};
		double[] tau = new double[]{0, 0, 0};
		int[] jpvt = new int[] {1, 2, 3};

		double[] VL = new double[M * N];
		double[] VR = new double[M * N];
		double[] WR = new double[N];
		double[] WI = new double[N];
		double[] S = new double[N];
		double[] U = new double[M * M];
		double[] VT = new double[N * N];
		double[] superb = new double[N - 1];
                int info=0;
                double[] work = new double[1000];
                int lwork = 3;
                int[] iwork = new int[24];
        //
        // Print the parameters
        //
		printMatrix("Matrix A", LUDecomposition.LAYOUT.RowMajor, A, M, K);
		printMatrix("Matrix B", LUDecomposition.LAYOUT.RowMajor, B, K, N);
        //
        // Set Option
        //
		int matrix_layout = LUDecomposition.LAYOUT.RowMajor;
		int trans = LUDecomposition.TRANSPOSE.NoTrans;
		int Uplo = LUDecomposition.UPLO.Upper;
        //
        //
		/* ---- LU ---- */
        //
		System.out.println("\n \n");
		System.out.println("###   LU   ###");
        //
		/* dcomp */
        //
		System.out.println();
		System.out.println("The LU decomposition A=PLU (L has unit diagnoal elements):");
		D = A.clone();
        int m = M;//m = the number of rows of the matrix a
        int n = N;//n = the number of columns of the matrix a
        int signum=1;
        double[] a = D;//a double precision, dimension(m,n)
        int lda = N;//lda = leading dimension of array a
        //ipiv = pivot indices  
        long[] p= new long[n];     
        signum=LUDecomposition.decomp( n, a, p, signum);
        printMatrix("  The LU matrix of A:",matrix_layout, a, n, n);
        for (int j=0; j<n; j++)
    	    System.out.printf("\t"+ string(p[j]));
        //
        /* dslve */
        //
        System.out.println();

	System.out.printf("The determinant of A: %.3f\n", LUDecomposition.det(n, a, signum));
	System.out.printf("The log-absolute value of the determinant of A: %.3f\n", LUDecomposition.lndet(n, a));

        D = B.clone();
        D2 = A.clone();
        a = D;
        signum=LUDecomposition.decomp( n, a, p, signum);

        //trans specifies the form of the system of equations
        n = M;//n = order of matrix a
        int nrhs = N;//nrhs = number of columns of matrix b 
        //a is computed by dgetrf
        //lda = n;//lda = leading dimension of a
        //ipiv = pivot indices from dgetrf
        double[] b = {12, 6, -4};//b double precision, dimension(n,nrhs)
        double[] x =new double[n];
        //int ldb = nrhs;//ldb = leading dimension of b
        LUDecomposition.slve( n, a, b, p, x);
        printMatrix("  The solution  matrix of BX=[12,6,4] by LU is:", matrix_layout, x, n, 1);
        
}


    /** Print the matrix X. */
    private static void printMatrix(String prompt, int layout, double[] X, int I, int J) {
    	System.out.println(prompt);
    	if (layout == LUDecomposition.LAYOUT.ColMajor) {
    		for (int i=0; i<I; i++) {
    			for (int j=0; j<J; j++)
    				System.out.print("\t" + string(X[j*I+i]));
    			System.out.println();
    		}
    	}
    	else if (layout == LUDecomposition.LAYOUT.RowMajor){
    		for (int i=0; i<I; i++) {
    			for (int j=0; j<J; j++)
    				System.out.print("\t" + string(X[i*J+j]));
    			System.out.println();
    		}
    	}
    	else{System.out.println("** Illegal layout setting");}
    }

    private static void printIntArray(String prompt, int[] X, int L) {
    	System.out.println(prompt);
    	for (int i=0; i<L; i++) {
    		System.out.print("\t" + string(X[i]));
    	}
    	System.out.println();
    }

    /** Shorter string for real number. */
    private static String string(double re) {
    	String s="";
    	if (re == (long)re)
    		s += (long)re;
    	else
    		s += re;
    	return s;
    }
}
