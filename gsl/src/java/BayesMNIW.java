package BayesMNIWConjugate;

import JAMAJniGsl.*;
import java.io.*;
import java.util.zip.GZIPInputStream;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;
import java.lang.*;


public class BayesMNIW implements java.io.Serializable{
    /* ------------------------
    * Class variables
    * ------------------------ */
    double v, vs;//scalar,v*, v*
    int n,q,p;//nsam means the number of random samples to generate
    double[] sigma, iden, zero;
    //double[][] beta,temp;
    Matrix X, Y, VL, Vrinv, Vs, Vsinv, Mub, Mus, Phi, Phis, SL, Beta, Iden;
    CholeskyDecomposition CholV,CholSi;
    final int matrix_layout = Matrix.LAYOUT.RowMajor;
    final int trans = Matrix.TRANSPOSE.Trans;
    final int notrans = Matrix.TRANSPOSE.NoTrans;
    /* ------------------------
    * Constructor
    * ------------------------ */
    public BayesMNIW( Matrix Yv, Matrix Xv, Matrix mubv, Matrix Vr, Matrix Phiv, double vv) {
    	/*Initiation*/
    	//nsam = nsamv;
    	Y = Yv;
    	X = Xv;
    	Mub = mubv;
    	Vrinv = Vr.inverse();
	Phi = Phiv;
	v = vv;
    	n = Y.getRowDimension();
	q = Y.getColumnDimension();
    	p = X.getColumnDimension();
    	sigma = new double[q * q];
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
	Vsinv = Vrinv.plus(X.transpose().times(X));
	//System.out.println(" OK");
	Vs = Vsinv.inverse();
	//System.out.println(" OK");
	Mus = Vs.times((X.transpose().times(Y)).plus(Vrinv.times(Mub)));
	//System.out.println(" OK");
	Phis = Phi.plus(Y.transpose().times(Y));
	//System.out.println(" OK");
	Phis = Phis.plus(Mub.transpose().times(Vrinv.times(Mub)));
	//System.out.println(" OK");
	Phis = Phis.minus(Mus.transpose().times(Vsinv.times(Mus)));
	//System.out.println(" OK");
	vs = v + (double)n;
        CholV = Vs.chol(); 
        VL = CholV.getL();

      	}
    /* ------------------------
    * Public Methods
    * ------------------------ */
    public Matrix getSigmaResult(){
	//Phis = Phis.inverse();
	double[] PhisArr = Phis.getRowPackedCopy();
	double[] result = new double[q * q];
	JniGslRng.inverwishart(vs,q,PhisArr,result);
	Matrix Result = new Matrix(result,q,q);
	//Result = Result.inverse();
    	return Result;
    }

    public Matrix getBetaResult(Matrix Msigma){
        CholSi = Msigma.chol(); 
        SL = CholSi.getL();	
	double[] sample = new double [p * q];
	JniGslRng.multinorm(p * q, zero,iden, sample);
	Matrix Sample = new Matrix(sample,p,q);
	Beta = VL.times(Sample);
	Beta = Beta.times(SL,notrans,trans);
	Beta = Beta.plus(Mus);
    	return Beta;
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
