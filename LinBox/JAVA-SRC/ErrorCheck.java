import java.io.*;
import java.net.URL;


/**
 * @author Brian Gold
 * @version 1.0
 * 
 * Date: 7/25/01
 * Last Updated: 1/4/02
 * File: ErrorCheck.java
 * 
 * This file contains only static methods which check user input for errors
 * For all methods, TRUE should be returned if there IS an error, and FALSE
 * should be returned if there IS NOT an error.
 */
public class ErrorCheck{
    

    
    /**
     * @param String email
     * @return boolean
     *
     * This method parses the string containing the email address
     * and determines if there is an error in the syntax.  If there is
     * an error, TRUE is returned.  If there is NOT an error, FALSE is returned
     * This method will only be called if the e-mail field is supposed to be
     * filled in.
     */
    public static boolean CheckEmailError(String email){

	// If the email field is blank, return an error
	if (email == null || email.equals("") ){
	    return true;
	}
	
	// If there are spaces in the email address, return an error
	else if (email.indexOf(' ') != (-1) ){
	    return true;
	}

	else{
	    // find the last occurrences of the symbols '@' and '.'
	    int indexAt = email.lastIndexOf('@');
	    int indexDot = email.lastIndexOf('.');

	    // If there is no '@' or no '.', return an error
	    if (indexAt == (-1) || indexDot == (-1) ){
		return true;
	    }

	    // If the last occurrence of '.' comes before the last occurrence
	    // of '@', return an error
	    else if (indexAt > indexDot){
		return true;
	    }
	}

	// else, the email is correctly formatted (Of course, it could
	// still be invalid.  If we are really concerned about this, 
	// we could try using the package JVerify available on the Internet
	return false;

    }  // end of method CheckEmailError

    
    /**
     * @param URLHandler muURL
     * @return boolean
     *
     * This method calls the method verifyURL() in file MatrixURL.java
     * to see if the URL is valid or not.  If it is not valid, TRUE is
     * returned, signifying an error.  If it is valid, this method returns
     * false, no error
     */
    public static boolean CheckURLError(URLHandler muURL){
	System.out.println("In ChecURLError");
	return (!muURL.verifyURL() );
    } // end of method CheckURLError
    

    /**
     * @param URLHandler muURL
     * @return boolean
     * 
     * This method calls the method getData() in the file MatrixURL.java
     * to perform a second check on the validity of the URL.  If it is NOT
     * VALID, TRUE is returned, signifying an error.  If it is valid, this
     * method returns false, no error
     */
    public static boolean CheckUrlDataError(URLHandler muURL)
	throws IOException{
	System.out.println("In CheckURLDataError");
	return (!muURL.getUrlData() );
    } // end of method CheckUrlDataError

    public static boolean CheckMatrixDataError(URLHandler muURL)
	throws IOException{
	System.out.println("In CheckMatrixDataError");
	return (!muURL.getMatrixData() );
    } // end of method CheckMatrixDataError


    /**
     * @param String stringCommands
     * @return boolean
     *
     * This method determines if there is anything in the text box
     * for user-entered commands in Gap.html besides the default
     * "RequirePackage("homology")".  If there is nothing more (or nothing
     * at all), TRUE is returned, signifying an error.  If there is other
     * text in the box, FALSE is returned
     */
    public static boolean CheckCommands(String stringCommands){
	if (stringCommands == null || stringCommands.equals("") ||
	    stringCommands.equals("RequirePackage(\"homology\");") )

	    return true;

	else
	    return false;
    }  // end of method CheckCommands






    /**
     * @param String execSearch
     * @return boolean
     *
     * This method determines if there is anything in the text box
     * for user-entered commands in Gap.html besides the default
     * "RequirePackage("homology")".  If there is nothing more (or nothing
     * at all), TRUE is returned, signifying an error.  If there is other
     * text in the box, FALSE is returned
     */
    public static boolean execSearch(String stringCommands){

	for(int i = 0; i < (stringCommands.length() - 5) ; i ++){

	    if (stringCommands.regionMatches(true, i, "Exec(", 0, 5) == true )
		return true;

	}

	    return false;
    }  // end of method execSearch





    /**
     * @param String procSearch
     * @return boolean
     *
     * This method determines if there is anything in the text box
     * for user-entered commands in Gap.html besides the default
     * "RequirePackage("homology")".  If there is nothing more (or nothing
     * at all), TRUE is returned, signifying an error.  If there is other
     * text in the box, FALSE is returned
     */
    public static boolean procSearch(String stringCommands){

	for(int i = 0; i < (stringCommands.length() - 8) ; i ++){

	    if (stringCommands.regionMatches(true, i, "Process(", 0, 8) == true )
		return true;

	}

	    return false;
    }  // end of method procSearch










    /**
     * @param String procSearch
     * @return boolean
     *
     * This method determines if there is anything in the text box
     * for user-entered commands in Gap.html besides the default
     * "RequirePackage("homology")".  If there is nothing more (or nothing
     * at all), TRUE is returned, signifying an error.  If there is other
     * text in the box, FALSE is returned
     */
    public static boolean commandCheck(String stringCommands){

	

	if (procSearch(stringCommands) || execSearch(stringCommands) )
		return true;

	else

	    return false;
    }  // end of method commandCheck





} // end of class ErrorCheck
