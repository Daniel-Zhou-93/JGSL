import JGSL.*;
import java.io.*;
import java.util.zip.GZIPInputStream;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;



public final class JGSLExamples {

	private JGSLExamples() {}

	public static void main(String[] args) {
        //
        // Prepare the matrices and other parameters
        //
        System.out.println("###   Exaples for JAMAJniGsl   ### \n \n");
        System.out.println("###   Parameter Preparation   ###");
        
        int M=3, N=3;
        double alpha = 2;
        int matrix_layout = Matrix.LAYOUT.RowMajor;
        long[] pivot;

        int trans = Matrix.TRANSPOSE.Trans;
    	int notrans = Matrix.TRANSPOSE.NoTrans;

        double[][] a = new double[][] {{12,-51,4},{6,167,-68},{-4,24,-41}};
        double[] a1 = new double[] {12, 6,-4};
        double[][] b = new double[][] {{4,-2,-6},{-2,10,9},{-6,9,14}};
        double[] b1 = new double[] {4,-2,-6};
        double[][] c;
        double[][] aa = new double[][] {{1,2,3},{2,4,4}};
        double[][] bb = new double[][] {{1,2},{3,1},{-1,0}};
        //
        //Construct Matrix
        //
        Matrix A = new Matrix(a);
        Matrix B = new Matrix(b);
        Matrix B1 = new Matrix(b1, N);
        Matrix A1 = new Matrix(a1, N);
        Matrix AA = new Matrix(aa);
        Matrix BB = new Matrix(bb);
	
	Matrix AT = new Matrix(A.getColumnDimension(), A.getRowDimension());
	Matrix AAT = new Matrix(AA.getColumnDimension(), AA.getRowDimension());
       	Matrix BBT = new Matrix(BB.getColumnDimension(), BB.getRowDimension());

        Matrix C = new Matrix(M,N);
        Matrix D = new Matrix(M,N);
        Matrix X = new Matrix(M,N);
        Matrix Q = new Matrix(M,N);
        Matrix R = new Matrix(M,N);
        Matrix L = new Matrix(M,N);
        Matrix U = new Matrix(M,N);
        Matrix V = new Matrix(M,N);
        Matrix S = new Matrix(M,N);
		//
        // Print the parameters
        //
        printMatrix("Matrix A", matrix_layout, A.getArray(), M, N);
        printMatrix("Matrix B", matrix_layout, B.getArray(), M, N);
        printMatrix("Matrix B1", matrix_layout, B1.getArray(), N, 1);
        printMatrix("Matrix A1", matrix_layout, A1.getArray(), N, 1);
        printMatrix("Matrix AA", matrix_layout, AA.getArray(), 2, 3);
	printMatrix("Matrix BB", matrix_layout, BB.getArray(), 3, 2); 

	System.out.println("M =" + M);
        System.out.println("N =" + N);
        System.out.println("alpha =" + alpha);
        //
        /* ---- Basic matrix algebra ---- */
        //
        System.out.println("\n\n###   Basic matrix algebra   ###");
        System.out.println();
        //
        //Matrix Addition
        //
        System.out.println("\n##  Addition: C = A + B  ##");
        C = A.plus(B);
        printMatrix("C = ", matrix_layout, C.getArray(), M, N);
        //
        //MINUS
        //
        System.out.println("\n##  Subtraction: C = A - B  ##");
        C = A.minus(B);
        printMatrix("C = ", matrix_layout, C.getArray(), M, N);
        //
        //TIMES
        //
        System.out.println("\n##  Multiplication: C = A * B  ##");
        C = A.times(B);
        printMatrix("C = ", matrix_layout, C.getArray(), M, N);
        //
        //TIMES with transpose option
        //
        System.out.println("\n##  Multiplication with transpose option: C = A' * B  ##");
        C = A.times(B,trans,notrans);
        printMatrix("C = ", matrix_layout, C.getArray(), M, N);
        //
        //TIMES with transpose option
        //
        System.out.println("\n##  Multiplication with transpose option: C = A * B  ##");
        C = A.times(B,notrans,notrans);
        printMatrix("C = ", matrix_layout, C.getArray(), M, N);
		//
        //TIMES with transpose option
        //
		printMatrix("Matrix AA", matrix_layout, AA.getArray(), 2, N);
        printMatrix("Matrix BB", matrix_layout, BB.getArray(), M, 2);
        System.out.println("\n##  Multiplication with transpose option: C = AA * BB  ##");
        C = AA.times(BB,notrans,notrans);
        printMatrix("C = ", matrix_layout, C.getArray(), C.getRowDimension(), C.getColumnDimension());
        //
        //TIMES with transpose option
        //
        System.out.println("\n##  Multiplication : C = AA' * BB'  ##");
        C = AA.times(BB, trans, trans);
        printMatrix("C = ", matrix_layout, C.getArray(), C.getRowDimension(), C.getColumnDimension());
        //
        //SCALAR
        //
        System.out.println("\n##  Scalar multiplication: C = alpha * A  ##");
        C = A.times(alpha);
        printMatrix("C = ", matrix_layout, C.getArray(), M, N);

	System.out.println("\n##  Scalar multiplication: C = alpha * A  ##");
        C = A.times(alpha);
        printMatrix("C = ", matrix_layout, C.getArray(), M, N);

	/* Transposition */
	System.out.println("\n##  Transpose AA' using gsl_transpose_memcpy function  ##");
        AAT = AA.TransposeMemcpy();
        printMatrix("AA' = ", matrix_layout, AAT.getArray(), AAT.getRowDimension(), AAT.getColumnDimension());

	System.out.println("\n##  Transpose square matrix A using gsl_transpose function  ##");
        AT = A.Transpose();
        printMatrix("AT = ", matrix_layout, AT.getArray(), AT.getRowDimension(), AT.getColumnDimension());

	}
    
	/* Print the matrix X */
    
    private static void printMatrix(String prompt, int layout, double[][] X, int I, int J) {
        System.out.println(prompt);
        if (layout == Matrix.LAYOUT.ColMajor) {
            for (int i=0; i<I; i++) {
                for (int j=0; j<J; j++)
                    System.out.print("\t" + string(X[j][i]));
                System.out.println();
            }
        }
        else if (layout == Matrix.LAYOUT.RowMajor){
            for (int i=0; i<I; i++) {
                for (int j=0; j<J; j++)
                    System.out.print("\t" + string(X[i][j]));
                System.out.println();
            }
        }
        else{System.out.println("** Illegal layout setting");}
    }

    /* Shorter string for real number */
    //
    private static String string(double re) {
        String s="";
        if (re == (long)re)
            s += (long)re;
        else
            s += re;
        return s;
    }


}
