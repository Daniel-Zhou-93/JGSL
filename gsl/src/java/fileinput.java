/*
 * fileinput.java
 *
 * Created on July 2, 2002, 8:09 PM
 */

package JGSL;

/**
 *
 * @author  sudipto banerjee
 * @version 
 */

import java.io.*;
import JGSL.*;
import java.util.StringTokenizer;

public class fileinput {

    /** Creates new fileinput */
    public fileinput() {
    }
    
    public static void readmatrix(File fin, Matrix A)
    {
     String temp;    
        try
        {
        FileReader fr = new FileReader(fin);
	BufferedReader br = new BufferedReader(fr);

        for (int i = 0; i < A.getRowDimension(); i++) {
            temp = br.readLine();
            StringTokenizer temp2 = new StringTokenizer(temp);//Comment: It's better to use "split" instead of Tokenizer
            //StringTokenizer temp2 = new StringTokenizer(temp,"	");// UPDATED ONE.
            for (int j = 0; j < A.getColumnDimension(); j++) {
		A.set(i,j,Double.parseDouble(temp2.nextToken()));
            }
        }
        }
        catch (IOException e) {
		System.out.println(e);
        }
    }
    
    public static void readvector(File fin, Matrix A)
    {
        String temp;     
        try
        {
        FileReader fr = new FileReader(fin);
	BufferedReader br = new BufferedReader(fr);

        for (int i = 0; i < A.getRowDimension(); i++) {
            temp = br.readLine();
            StringTokenizer temp2 = new StringTokenizer(temp);
            A.set(i,0,Double.parseDouble(temp2.nextToken()));            
        }
        }
        catch (IOException e) {
            System.out.println(e);
        }
    }
    
    public static void readArray(File fin, double [][] A, int nrows, int ncols)
    {
        String temp;     
        try
        {
        FileReader fr = new FileReader(fin);
	BufferedReader br = new BufferedReader(fr);

        for (int i = 0; i < nrows; i++) {
            temp = br.readLine();
            StringTokenizer temp2 = new StringTokenizer(temp);
            for (int j = 0; j < ncols; j++) {
		A[i][j] = Double.parseDouble(temp2.nextToken());
            }
        }
        }
        catch (IOException e) {
		System.out.println(e);
        }
    }
    
    public static void readArray(File fin, int [][] A, int nrows, int ncols)
    {
        String temp;     
        try
        {
        FileReader fr = new FileReader(fin);
	BufferedReader br = new BufferedReader(fr);

        for (int i = 0; i < nrows; i++) {
            temp = br.readLine();
            StringTokenizer temp2 = new StringTokenizer(temp);
            for (int j = 0; j < ncols; j++) {
		A[i][j] = Integer.parseInt(temp2.nextToken());
            }
        }
        }
        catch (IOException e) {
		System.out.println(e);
        }
    }
    
    public static void readArray(File fin, double [] A)
    {
        String temp;     
        try
        {
        FileReader fr = new FileReader(fin);
	BufferedReader br = new BufferedReader(fr);

        for (int i = 0; i < A.length; i++) {
            temp = br.readLine();
            StringTokenizer temp2 = new StringTokenizer(temp);
		A[i] = Double.parseDouble(temp2.nextToken());
            
        }
        }
        catch (IOException e) {
		System.out.println(e);
        }
    }
    
    public static void readArray(File fin, int [] A)
    {
        String temp;     
        try
        {
        FileReader fr = new FileReader(fin);
	BufferedReader br = new BufferedReader(fr);

        for (int i = 0; i < A.length; i++) {
            temp = br.readLine();
            StringTokenizer temp2 = new StringTokenizer(temp);
		A[i] = Integer.parseInt(temp2.nextToken());
            
        }
        }
        catch (IOException e) {
		System.out.println(e);
        }
    }
    
    
}
