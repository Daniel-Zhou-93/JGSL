import JGSL.*;
import java.io.*;
import java.util.zip.GZIPInputStream;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;


public final class BayesTest {
	private BayesTest() {}
	public static void main(String[] args) {
        // Test NIG
        // Generate one Normal-Inverse Gamma random variable.
        int matrix_layout = Matrix.LAYOUT.RowMajor;
	JniGslRng.seed();
        double[] mu_sigma2 = new double[2];
        double mu = 2.0;
        double lambda = 1.5;
        double a = 3;
        double b = 5;
        mu_sigma2 = BayesNIG.rnorminvergamma(mu, lambda, a, b);
        System.out.printf("One-sample Normal Inverse Gamma: mu %.2f, sigma2 %.2f\n", mu_sigma2[0], mu_sigma2[1]);

        // Generate three Normal-Inverse Gamma random variables in one collection.
        double[][] mu_sigma2_x3 = new double[3][2];
        mu_sigma2_x3 = BayesNIG.rmvnorminvergamma(3, mu, lambda, a, b);
        System.out.printf("Three-sample Normal Inverse Gamma (True mu: %f):\n", mu);
        System.out.printf("Sample 1: mu %.2f, sigma2 %.2f\n", mu_sigma2_x3[0][0], mu_sigma2_x3[0][1]);
        System.out.printf("Sample 2: mu %.2f, sigma2 %.2f\n", mu_sigma2_x3[1][0], mu_sigma2_x3[1][1]);
        System.out.printf("Sample 3: mu %.2f, sigma2 %.2f\n", mu_sigma2_x3[2][0], mu_sigma2_x3[2][1]);

        // Test NIW
        int q = 4;
	double[][] PhivPrior = {{10.3838, 1.06197, -2.65679, -2.84745},{1.06197, 8.84605, 0.122342, 0.746835},{-2.65679,0.122342, 3.77193, -2.28483},{-2.84745, 0.746835, -2.28483, 10.571}};
	Matrix Phiv = 	new Matrix(PhivPrior);//phi q*q
        double n = 7;
        Matrix Sigma = new Matrix(q, q);
        Sigma = BayesNIW.getSigmaResult(q, n, Phiv);
        double[] muv = {1.0, 2.0, 3.0, -1.0};
        double[] mu_out = new double[q];

        System.out.println("Phi:");
        printMatrix(matrix_layout, Phiv.getArray(), Phiv.getRowDimension(), Phiv.getColumnDimension());
        System.out.println("Sigma:");
        printMatrix(matrix_layout, Sigma.getArray(), Sigma.getRowDimension(), Sigma.getColumnDimension());

        // subtract Sigma from its symmetric view
        System.out.println("Internal Debug: Get difference with symmetric matrix.");
        double[][] symmSigma = Sigma.getArrayCopy();
        for (int i = 0; i < q; i++) {
            for (int j = i + 1; j < q; j++) symmSigma[j][i] = symmSigma[i][j];
        }
        
        Matrix symmSigmaMat = new Matrix(symmSigma);
        Matrix deltaSigma = Sigma.minus(symmSigmaMat);
        printMatrix(matrix_layout, deltaSigma.getArray(), Sigma.getRowDimension(), Sigma.getColumnDimension());
        System.out.println("Symmetric difference is small. Set Sigma to its symmetry.");
        
        Sigma = symmSigmaMat;
        mu_out = BayesNIW.getMuResult(q, muv, lambda, Sigma);
        printDoubleArray("Actual Mu:", muv, q);
        printDoubleArray("Sampled Mu:", mu_out, q);
        
        // get a pds matrix to simulate sampling from a multivariate normal.
        
        Matrix V = new Matrix((int)n, (int)n);
        String vfilename = "../data/V.txt";
        File vfile = new File(vfilename);
        fileinput vfi = new fileinput();
        vfi.readmatrix(vfile, V);
        
        double[][] Mu_array = {{1,0,1,2},{2,0,2,3},{3,0,3,5},{5,2,5,9},{8,17,7,-8},{13,20,4,9},{21, 8, 3, -4}};
        Matrix Mu_s = new Matrix(Mu_array);
        Matrix Mu_out = new Matrix((int)n, q);
        Mu_out = BayesNIW.getMatrixMuResult((int)n, q, Mu_s, Sigma, V);
        System.out.println("Sample 4 x 7 matrix Mu, drawn from a multivariate normal.");
        printMatrix(matrix_layout, Mu_out.getArray(), Mu_out.getRowDimension(), Mu_out.getColumnDimension());
        
	JniGslRng.free();
        
    }

    //
    /* Print the matrix X */
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
    private static void printMatrix(int layout, double[][] X, int I, int J) {
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
    /* Print X but it's a double array */
    //
    private static void printDoubleArray(String prompt, double[] X, int L) {
    	System.out.println(prompt);
    	for (int i=0; i<L; i++) {
    		System.out.printf("\t%.3f", X[i]);
    	}
    	System.out.println();
    }

    //
    /* Print X, but it's a generic array. */
    // TODO: under development
    //
    private static <T> void printArray(String prompt, T[] X, int L) {
    	System.out.println(prompt);
    	for (int i=0; i<L; i++) {
    		System.out.printf("\t%.3f", X[i]);
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
}
