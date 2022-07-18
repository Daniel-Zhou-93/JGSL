package JGSL;

import JGSL.*;
import java.io.*;
import java.util.zip.GZIPInputStream;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;
import java.lang.*;

/*
Constructor for normal inverse gamma random variables.
*/

public class BayesNIG {
    /* ------------------------
    * Class variables
    * ------------------------ */
    double[] z;
    double[] sigma2;
    final int matrix_layout = Matrix.LAYOUT.RowMajor;
    final int trans = Matrix.TRANSPOSE.Trans;
    final int notrans = Matrix.TRANSPOSE.NoTrans;
    /* ------------------------
    * Constructor
    * ------------------------ */
    // NOTE: Don't take Yv and Xv. Keep others. Sample inverse gamma 1/sigma^2 from gamma (a,b)
    public BayesNIG() {}
    //
    /* Random number generator; generates a single Normal-Inverse Gamma pair. */
    //
    public static double[] rnorminvergamma(double mu, double lambda, double a, double b) {
        double[] mu_sigma2 = new double[2];
        mu_sigma2[1] = JniGslRng.invergamma(a, b);
        mu_sigma2[0] = mu + JniGslRng.uninorm(Math.sqrt(mu_sigma2[1]/lambda));
        return mu_sigma2;
    }

    //
    /* Random number generator; generates multiple NIG variables */
    //
    public static double[][] rmvnorminvergamma(int nsam, double mu, double lambda, double a, double b) {
        double[][] mu_sigma2s = new double[nsam][2];
     	for(int i=0; i<nsam; i++){
            // independently sampling with a different inverse gamma each, but same mu
            mu_sigma2s[i] = rnorminvergamma(mu, lambda, a, b);
     	}
        return mu_sigma2s;
//    	Matrix Z = new Matrix(z,p);
//    	mu = mu.plus(L.times(sigma).times(Z));
//	    return mu.getRowPackedCopy();
    }
}
