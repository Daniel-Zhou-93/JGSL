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
    double v, vs;//scalar,v*, v*
//    int n,q,p;//nsam means the number of random samples to generate
    double[] iden, zero; //mu, sigma, iden, zero;
    //double[][] beta,temp;
    Matrix VL, Vs, Vrinv, Mub, Mus, Phi, SL, Mu_out, Iden;
    CholeskyDecomposition CholV,CholSi;
    final int matrix_layout = Matrix.LAYOUT.RowMajor;
    final int trans = Matrix.TRANSPOSE.Trans;
    final int notrans = Matrix.TRANSPOSE.NoTrans;
    /* ------------------------
    * Constructor
    * ------------------------ */
    public BayesNIW() {}

    /* ------------------------
    * Another Constructor to store variables, for easier management of variables without having to allocate
    * ------------------------ */
    public BayesNIW(int p, int q, Matrix mubv, Matrix Vr, Matrix Phiv, double vv) {
    	/*Initiation*/
        /*p is the side length of row variance matrix Vr, q is the side length of Phiv */
    	//nsam = nsamv;
    	Mub = mubv;
        CholV = new CholeskyDecomposition(Vr); 
        VL = CholV.getL();
        Vrinv = Vr.copy();
        double[] VrArray = Vrinv.getRowPackedCopy();
        Vrinv = new Matrix(VrArray, p, p);
	Phi = Phiv;
	v = vv;
//    	n = Y.getRowDimension();
//	q = Y.getColumnDimension();
//    	p = X.getColumnDimension();
//    	sigma = new double[q * q];
    	//sigma2 = new double[nsam];
    	//sigma2recip = new double[nsam];
    	//beta = new double[nsam][p];
	Iden = Matrix.identity(p * q, p * q);
	iden = Iden.getRowPackedCopy();
	zero = new double[p * q];
	//double[][] XXX = {{1,2,3},{1,2,3}};
	//Matrix XXXX = new Matrix(XXX);
	//Matrix XXXXX = XXXX.times(XXXX,trans,notrans);
	//System.out.println(" OK");
    	/*Calculation*/
	//Vsinv = Vrinv.plus(X.times(X,trans,notrans));
//	Vsinv = Vrinv.plus(X.transpose().times(X));
	//System.out.println(" OK");
//	Vs = Vsinv.inverse();
	//System.out.println(" OK");
//	Mus = Vs.times((X.transpose().times(Y)).plus(Vrinv.times(Mub)));
	//System.out.println(" OK");
//	Phis = Phi.plus(Y.transpose().times(Y));
	//System.out.println(" OK");
//	Phis = Phis.plus(Mub.transpose().times(Vrinv.times(Mub)));
	//System.out.println(" OK");
//	Phis = Phis.minus(Mus.transpose().times(Vsinv.times(Mus)));
	//System.out.println(" OK");
//	vs = v + (double)n;
//        CholV = Vs.chol(); 
//        VL = CholV.getL();

      	}

    /* ------------------------
    * Public Methods
    * ------------------------ */
    public static Matrix getSigmaResult(int q, double n, Matrix Phi){
	//Phis = Phis.inverse();
        CholeskyDecomposition CholPhi = new CholeskyDecomposition(Phi); 
	double[] PhiArr = CholPhi.getL().getRowPackedCopy();
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
    //* TODO: write MNIW function to sample a MNIW matrix */
    //
    public static Matrix getMuResult(int p, int q, Matrix mumat, double lambda, Matrix Msigma, Matrix MV){
        //TODO: sample from matrix normal distribution.
        CholeskyDecomposition CholSi = new CholeskyDecomposition(Msigma.times(1.0/lambda)); 
        Matrix mu_out = mumat;//new Matrix();
        return mu_out;
    }
    // TODO: copy this function as a static function and then produce 
    // TODO: initialize values with constructor. You can also produce a non-constructor variant, so you don't have to keep allocating.
    //* Sample a matrix normal variable given row covariance matrix V and column covariance matrix MSigma. */
    public static Matrix getMatrixMuResult(int p, int q, Matrix Mus, Matrix Msigma, Matrix V){
        int trans = Matrix.TRANSPOSE.Trans;
        int notrans = Matrix.TRANSPOSE.NoTrans;
        CholeskyDecomposition CholSi = new CholeskyDecomposition(Msigma); 
        Matrix SL = CholSi.getL();	
	double[] sample = new double[p * q];
        Matrix Iden = Matrix.identity(p * q, p * q);
        double[] iden = Iden.getRowPackedCopy();
        double[] zero = new double[p * q];
	JniGslRng.multinorm(p * q, zero,iden, sample);
	Matrix Sample = new Matrix(sample,p,q);
        CholeskyDecomposition CholV = new CholeskyDecomposition(V);
        Matrix VL = CholV.getL();
	Matrix Mu_out = VL.times(Sample);
	Mu_out = Mu_out.times(SL,notrans,trans);
	Mu_out = Mu_out.plus(Mus);
    	return Mu_out;
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
