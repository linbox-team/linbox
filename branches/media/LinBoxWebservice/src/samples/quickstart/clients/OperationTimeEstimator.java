// OperationTimeEstimator.java
// Matthew Fendt and David Saunders, Summer '08
// Provides simple time estimates for the duration of certain matrix operations



package samples.quickstart.clients;

import java.lang.*;

public class OperationTimeEstimator
{
    
    /*
      Argument meanings:
      
      m = number of rows.
      n = number of cols.
      s = number of entries specified in input data file (any other matrix 
          entries are implicitly zero).
      d = number of bits in largest entry, approximately 
          max_v(1 + log_2(abs(v))).
      r = rank, 0 <= r <= min(m, n).  Generally unknown to user, but 
          occasionally known a priori.
      

      File formats:
       SMS format example:
       4 6 M
       1 1 3
       1 3 7
       2 2 9
       3 1 27
       4 3 -100
       0 0 0
       m = 4, n = 6, s = 5, d = 1 + log(100), r = 4.
       s is from number of lines strictly between '4 6 M' and '0 0 0'.
       d from 3 7 9 27 -100
       
       Other basic sparse format:
       4 6 5
       1 1 3
       1 3 7
       2 2 9
       3 1 27
       4 6 -100
       parameters as before, s is from third number on line 1.
       
       Dense format:
       4 6
       3  0 7 0 0  0
       0  9 0 0 0  0
       27 0 0 0 0  0
       0  0 0 0 0 -100
       m = 4, n = 6, s = 24, d = 1 + lg(100),
       s is from m*n.
    */





    // ************************************************************************
    // ****** Machine specific costs.

    // true iff machine cost factors have been obtained.
    static boolean GoodMachineCosts = false;
    static double C; // cost of wordsize modular axpy.
    static double D; // cost per entry of wordsize modular dot product.
    // Thus D*n is cost of length n vector dot product, wordsize modular.

    // average cost of wordsize modular axpy, when implicit within blas.
    static double E;

    int MAX_LENGTH = 100;
    // ************************************************************************

    static void initMachineCosts (/* url Service */) 
	// to become a service call getting these values.
    {   C = 1/25.42e6; D = 1/20.46e6; E = 0.5e-8; GoodMachineCosts = true;   }
    // ************************************************************************


    boolean isDense (int m, int n, double s) 
	// a kludge for now.
    {   return m*n < 3*s;   }



    // *********************************************************************** 

    // cost estimate in seconds to compute minimum polynomial of mxn integer 
    // matrix of rank r and having s nonzero entries of size bounded by d bits.
    double minpolyCost (int m, int n, double s, int d, int r)    
    {
	if ( ! GoodMachineCosts ) initMachineCosts();
	if (m > n) 
	    {
		// Swamp m and n
		int temp = m;
		m = n;
		n = temp;
	    }

	if (r == 0) r = m;
	double F = 1; // minpoly fudge factor.

	return isDense(m,n,s) ? Math.ceil(F*E*d*m*n*r) : 
	    Math.ceil(F*d*(C*s + D*n)*r);
    }

    // ************************************************************************
    double detCost (int m, int n, double s, int d)       
    {
	int r = 0;

	if ( ! GoodMachineCosts ) initMachineCosts();

	if (m > n)
	    {
		// Swap m and n
		int temp = m;
		m = n;
		n = temp;
	    }

	if (r == 0) r = m;
	double F = 1; // det fudge factor.
	
	return isDense(m,n,s) ? Math.ceil(F*d*E*m*n*r) : 
	    Math.ceil(F*d*minpolyCost(m, n, s, d, m));
    }

    // ************************************************************************
    double rankCost (int m, int n, double s, int d)
    { 
	int r = 0;

	if ( ! GoodMachineCosts ) initMachineCosts();

	if (m > n)
	    {
		// Swap m and n
		int temp = m;
		m = n;
		n = temp;
	    }

	if (r == 0) r = m;
	double F = 2; // rank fudge factor.
	return isDense(m,n,s) ? Math.ceil(F*E*m*n*r) : 
	    Math.ceil(F*minpolyCost(m, n, s, 1, m));
    }

    // ************************************************************************
    double valenceCost (int m, int n, double s, int d)
    { 
	double F = 1; // valence fudge factor.
	return Math.ceil(F*minpolyCost(m,n,s,d, 0)); 
    }


    // ************************************************************************


    double traceCost (int m, int n, double s, int d)
    {	// Trace is very fast.  Perhaps just return 1 second, ignoring input.
	double F = 0.5; // trace fudge factor.
	
	if (m > n) 
	    {
		// Swamp m and n
		int temp = m;
		m = n;
		n = temp;
	    }
	
	return Math.ceil(F*C*d*m);
    }

    // ************************************************************************

    double smithformCost (int m, int n, double s, int d)
    {
	int r = 0;

	double F = 1.0/6; // smithform fudge factor.

	if (m > n) 
	    {
		// Swamp m and n
		int temp = m;
		m = n;
		n = temp;
	    }

	if (r == 0) r = m;
	
	// cubic when m = n = r.  When structure is rich expect m^3.5 or even m^4.
	return Math.ceil(F*E*d*m*n*r);  // future tuning...
    }

 // ************************************************************************


    double solveCost (int m, int n, double s, int d)
    {
	int r = 0;
	double F = 1; // solve fudge factor.
	double BBfactor = 6;

	if (m > n) 
	    {
		// Swamp m and n
		int temp = m;
		m = n;
		n = temp;
	    }

	if (r == 0) r = m;
	
	return isDense(m,n,s) ? Math.ceil(F*d*E*m*n*r) : 
	    Math.ceil(BBfactor*F*d*minpolyCost(m, n, s, d, m));
    }

 // ************************************************************************





    String prettyIter (double t, int k)
    {
	int m = 60, h = 60*m, d = 24*h;
	StringBuffer foo = new StringBuffer(MAX_LENGTH);
	if (t < 1) foo.append("1s");
	else if (k > 1)
	    {
		if (t < m)
		    { foo.append(Math.floor(t)); foo.append("s");}

		else if (t < h) 
		    {foo.append(Math.floor(t/m)); foo.append("m "); 
		    foo.append(prettyIter(t - m*Math.floor(t/m),1)); }

		else if (t < d) 
		    {foo.append(Math.floor(t/h)); foo.append("h ");
		    foo.append(prettyIter(t - h*Math.floor(t/h),1));}

		else 
		    {foo.append(Math.floor(t/d)); foo.append("d "); 
		    foo.append(prettyIter(t - d*Math.floor(t/d),1));}
	    }
	else 
	    {
		if (t < m)
		    { foo.append(Math.ceil(t)); foo.append("s");}

		else if (t < h) 
		    {foo.append(Math.ceil(t/m)); foo.append("m");}

		else if (t < d) 
		    {foo.append(Math.ceil(t/h)); foo.append("h");}

		else 
		    {foo.append(Math.ceil(t/d)); foo.append("d");}
	    }

	return foo.toString();
    }


    // ************************************************************************

    String pretty(double t) { return prettyIter (t, 2); }

    // ************************************************************************

    public OperationTimeEstimator()
    {
	initMachineCosts();
    }


    // ************************************************************************


    public String estimateTime(String operation, String matrix, 
			       boolean detailedEstimate)
    {
	try {
	// The estimated time
	double estimate;

	// Delineators for the matrix string
	String[] temp = matrix.split("%0A|%0D|\\+");

	// The rows and columns of the matrix
	double m = -1;
	double n = -1;

	boolean pastM = false;
	boolean firstSet = true;
	boolean dense = false;


	int count = 1;

	// The largest entry
	double d = -1;

	// The number of non zero entries
	int s = 0;
	
	for (int j = 0; j < temp.length; j++)	
	    {	   
		//System.out.println(temp[j]);

		// Set m
		if (m == -1 && !temp[j].equals(""))
		    {
			m = Double.parseDouble(temp[j]);
			System.out.println("M is " + m);
		    }
		
		// Set n
		else if (n == -1 && m != -1 && !temp[j].equals(""))
		    {
			n = Double.parseDouble(temp[j]);
			System.out.println("N is " + n);
		    }
		
		// See if there is an 'M' or 'S' in the matrix file
		else if (pastM == false && m != -1 && n != -1)

		    {
			// Sparse format
			if (temp[j].equalsIgnoreCase("m") ||
			  temp[j].equalsIgnoreCase("s"))
			    {
				pastM = true;
				dense = false;
			    }

			// Dense format
			else
			    {
				pastM = true;
				dense = true;
			    }
		    }

		// Find the largest value
		else if (pastM == true && !temp[j].equals(""))
		    {
			// Since the file is in <i, j, v> pairs, we
			// only need to look at every third number to
			// determine the largest value
			if (count < 3 && dense == false)
			    count ++;
			
			// Temporarily let the first value in the
			// matrix be the largest
			else if (firstSet == true)
			    {

				firstSet = false;

				d = Math.abs(Double.
					     parseDouble(temp[j]));
				count = 1;
				
				s++;
			    }
			
			// If this value is bigger than the largest,
			// it is the new largest
			else if (firstSet == false)
			    {
				double tempD = 
				    Double.parseDouble(temp[j]);

				if (Math.abs(tempD) > d)
				    d = Math.abs(tempD);
				
				count = 1;
				
				s++;
			    }
		    }			
	    }
	    

	// Take log base 2 of d
	Double tempDouble = Math.log(d) / Math.log(2);
	d = tempDouble.intValue() + 1;


	int d1 = (int)d;
	int m1 = (int)m;
	int n1 = (int)n;

	System.out.println("D is " + d1);
	System.out.println("S is " + s);


	if (operation.equalsIgnoreCase("rank"))
	    estimate = rankCost(m1,n1,s,d1);
	
	else if (operation.equalsIgnoreCase("determinant"))
	    estimate = detCost(m1,n1,s,d1);

	else if (operation.equalsIgnoreCase("valence"))
	    estimate = valenceCost(m1,n1,s,d1);

	else if (operation.equalsIgnoreCase("trace"))
	    estimate = traceCost(m1,n1,s,d1);

	else if (operation.equalsIgnoreCase("smithNormalForm"))
	    estimate = smithformCost(m1,n1,s,d1);

	else
	    return "No time estimate available for " + operation;

	// The answer as a readable string	
	if (detailedEstimate == true)
		return pretty(estimate);

	// The answer in seconds
	else 
	    return Double.toString(estimate);
	}
	
	catch (Exception e)
	    {
		return "-1";
	    }
    }
}


