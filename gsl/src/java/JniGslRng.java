package JGSL;

/**gslRng.java*/

public class JniGslRng {
 private JniGslRng() {}

 static {
     
    /* load library (which will contain wrapper for Rmath function.)*/
    System.loadLibrary("gsl_random");
 
 }

	/* Seed */

public static native void	seed();


	/* Normal Distribution */

public static native double	uninorm(double sigma);
public static native void	multinorm(int k, double[] mu, double[] l, double[] result);

	/* Gamma Distribution */

public static native double	gamma(double a, double b);
public static native double	invergamma(double a, double b);


	/* wishart Distribution */

public static native void	wishart(double n, int p, double[] l, double[] reulst);
public static native void	inverwishart(double n, int p, double[] l, double[] result);


	/* free */

public static native void	free();

    
}







