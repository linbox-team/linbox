package LinBoxServerInterface;

import java.io.File;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.IOException;
import java.util.Random;

/** A class with only static methods that will generate matrix files in sms
  * format with given specifications.
  */
public class MatrixGenerator implements Locations {
	static private File getFile( int id ) {
		return new File( matrixDir, String.valueOf( id ) + ".sms" );
	}

	static private PrintWriter makeWriter( File f, int rows, int cols )
		throws IOException
	{	PrintWriter out = new PrintWriter (
			new BufferedWriter( new FileWriter(f) )
			                          );
		out.println( rows + " " + cols + " M" );
		return out;
	}

	static private void writeFooter( PrintWriter out )
		throws IOException
	{
		out.println( "0 0 0" );
		out.close();
	}

	/** Generate a new random matrix.
	  * @param id The id of the matrix, used to generate the file.
	  * @param rows The number of rows in the matrix
	  * @param cols The number of columns in the matrix
	  * @param sparsity The sparsity of the matrix.  This is defined as the
	  *                 probability that a given cell is nonzero, NOT the
	  *                 percentage of cells that are nonzero.  Should be a
	  *                 value between zero and one, inclusive.
	  * @return The location of the generated matrix file.
	  */
	static public File newRandom( int id, int rows, int cols, 
	                              float sparsity )
	{
	    File toRet = getFile( id );
	    try {
		PrintWriter out = makeWriter( toRet, rows, cols );
		Random generator = new Random();
		for( int i = 1; i <= rows; ++i )
		    for( int j = 1; j <= rows; ++j ) {
			if( generator.nextFloat() < sparsity )
			    out.println( i + " " + j + " " +
					 generator.nextInt() );
		    }
		writeFooter( out );
	    }
	    catch( IOException e ) {
	    	throw new RuntimeException( e );
	    }
	    return toRet;
	}

	/** Generate a new zero matrix.
	  * @param id The id of the matrix, used to generate the file.
	  * @param rows The number of rows in the matrix
	  * @param cols The number of columns in the matrix
	  * @return The location of the generated matrix file.
	  */
	static public File newZero( int id, int rows, int cols ) {
	    File toRet = getFile( id );
	    try {
		PrintWriter out = makeWriter( toRet, rows, cols );
		writeFooter( out );
	    }
	    catch( IOException e ) {
	    	throw new RuntimeException( e );
	    }
	    return toRet;
	}

	/** Generate a new identity matrix.
	  * @param id The id of the matrix, used to generate the file.
	  * @param rows The number of rows and columns in the matrix.  (Identity
	  *             matrices must be square.)
	  * @return The location of the generated matrix file.
	  */
	static public File newIdentity( int id, int rows ) {
	    File toRet = getFile( id );
	    try {
		PrintWriter out = makeWriter( toRet, rows, rows );
		for( int i = 1; i <= rows; ++i )
			out.println( i + " " + i + " " + "1" );
		writeFooter( out );
	    }
	    catch( IOException e ) {
	    	throw new RuntimeException( e );
	    }
	    return toRet;
	}
}
