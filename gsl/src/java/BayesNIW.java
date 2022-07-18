package JGSL;

import JGSL.*;
import java.io.*;
import java.util.zip.GZIPInputStream;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;
import java.lang.*;


public class BayesNIW implements java.io.Serializable{
    /* ------------------------
    * Class variables
    * ------------------------ */
//    double v, vs;//scalar,v*, v*
//    int n,q,p;//nsam means the number of random samples to generate
//    double[] mu, sigma, iden, zero;
    //double[][] beta,temp;
    // TODO: write single normal inverse wishart, not MNIW.
//    Matrix X, Y, VL, Vrinv, Vs, Vsinv, Mub, Mus, Phi, Phis, SL, Beta, Iden;
//    CholeskyDecomposition CholV,CholSi;
    final int matrix_layout = Matrix.LAYOUT.RowMajor;
    final int trans = Matrix.TRANSPOSE.Trans;
    final int notrans = Matrix.TRANSPOSE.NoTrans;
    /* ------------------------
    * Constructor
    * ------------------------ */
    public BayesNIW() {}

    /* ------------------------
    * Public Methods
    * ------------------------ */
    public static Matrix getSigmaResult(int q, double n, Matrix Phi){
	//Phis = Phis.inverse();
	double[] PhiArr = Phi.getRowPackedCopy();
	double[] result = new double[q * q];
	JniGslRng.inverwishart(n,q,PhiArr,result);
	Matrix Result = new Matrix(result,q,q);
	//Result = Result.inverse();
    	return Result;
    }

    public static double[] getMuResult(int q, double[] muv, double lambda, Matrix Msigma){
        CholeskyDecomposition CholSi = new CholeskyDecomposition(Msigma.times(1.0/lambda)); 
        Matrix SL = CholSi.getL();	
	double[] mu_out = new double[q];
	JniGslRng.multinorm(q, muv, SL.getRowPackedCopy(), mu_out);
        return mu_out;
    }

    //
    //* Print the matrix X */
    //
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
  	//
	/* Print the array X */
  	//
	private static void printIntArray(String prompt, int[] X, int L) {
  		System.out.println(prompt);
  		for (int i=0; i<L; i++) {
   		System.out.print("\t" + string(X[i]));
 	}
 	System.out.println();
	}

  	//
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

	private static final long serialVersionUID = 1;

}
