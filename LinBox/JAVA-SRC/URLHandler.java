import java.util.*;
import java.util.zip.*;
import java.net.*;
import java.io.*;
import java.text.*;


/**
 * @author Brian Gold
 * @version 1.0
 * 
 * Date: 7/25/01
 * Last Updated: 1/4/02
 * File: URLHandler.java
 * 
 * This file contains a class definition that takes a string and 
 * makes a URL out of it.  It also contains a String for data, the data
 * being all the information on the URL.
 */
public class URLHandler{

    // Data members of the URLHandler class
    public URL url;       // holds the URL of the site
    public String data;   // contains all the data of the site in the
                          // data member URL url

    /**
     * @param String stringURL
     *
     * This method is the constructor for the URLHandler class
     * which converts the String containing a url into a URL object
     * and initializes the data to ""
     */
    public URLHandler(String stringURL){

	// Try to extract URL from string
	try{
	    url = new URL(stringURL);
	}

	// If the URL is invalidly formed, set it to NULL
	catch(MalformedURLException murle){
	    url = null;
	}

	// Initialize the data
	data = "";

    } // End of constructor URLHandler


    /**
     * @param none
     * @return boolean
     *
     * This method returns true if the URL is correctly formed and valid;
     * otherwise, it returns false
     */
    public boolean verifyURL(){

	// if the URL is invalidly formed, return false
	if (url == null)
	    return false;

	// else, URL is correctly formed,
	// so now check the URL to see if it does exist
	try{
	    return doesURLExist();
	}
	catch(IOException ioe){
	    return false;
	}
    } // end of method verifyURL


    /**
     * @param none
     * @return boolean
     *
     * This method is called by verifyURL and actually checks to see
     * if the URL exists, since if it has gotten this far it is correctly
     * formed.  It returns true if the URL is valid, otherwise it returns false
     */
    public boolean doesURLExist() throws IOException{
	
	// Make a connection to the URL and get the length of the URL's data
	URLConnection siteCon = url.openConnection();
	int len = siteCon.getContentLength();

	// If the length is less than or equal to 0, the page is an invalid
	// page, so return false
	if (len <= 0)
	  return false;

	// else, the page is valid, so return true
	else
	  return true;

    } // End of method doesURLExist


    /**
     * @param none
     * @return boolean
     *
     * This method stores the data from the URL into the data member String
     * 'data'.  It returns true if there is data, and false if there is not
     */
    public boolean getUrlData() throws IOException{
        try{
            // Set up the Streams to read in the data from the URL
            InputStream uin = url.openStream();
            BufferedReader in = new BufferedReader(new InputStreamReader(uin));
            StringBuffer sb = new StringBuffer(2);
            String line;  // to hold each line read in, a temporary string

            // Read all the data in
            while (( line = in.readLine()) != null)
                sb.append(line).append("\n");
            
            data = sb.toString();
        }

	// If the URL is invalid, display an error (though it should be valid
        // if it reaches this point)
	catch(FileNotFoundException fnfe){
	    //System.out.println("Filenotfound in URLExtractor.java");
	    return false;
	}

	// If there is no data on the URL's page, return false
	if (data == null || data == ""){
	    //System.out.println("New false checker");
	    return false;
	}
	else
	    return true;
	
    } // End of method getUrlData


    /**
     * @param none
     * @return boolean
     *
     * This method stores the data from the URL into the file named by the
     * variable 'data'.  It returns true if there is data, and false if
     * there is not
     */
     public boolean getMatrixData() throws IOException{
	try{
	    // Set up the Streams to read in the data from the URL
	    InputStream uin = url.openStream();
	    BufferedReader in = new BufferedReader(new InputStreamReader(uin));
            
            File URLFILE = new File( url.getFile() );
            data = URLFILE.getName();
	    if (data.endsWith(".gz")==true){
	      try{
		GZIPInputStream zip = new GZIPInputStream(uin);
		FileOutputStream out = new FileOutputStream(data + ".txt");
		byte[] buffer = new byte[256];
		while (true) {
		  int bytesRead = zip.read(buffer);
		  if (bytesRead == -1) break;
		  out.write(buffer, 0, bytesRead);
		  //System.out.println(buffer);
		}
		out.close();
		FileInputStream fin = new FileInputStream(data + ".txt");
		FileOutputStream fout = new FileOutputStream(data);
		GZIPOutputStream gout = new GZIPOutputStream(fout);
		buffer = new byte[256];
		while(true){
		  int bytesRead2 = fin.read(buffer);
		  if (bytesRead2 == -1) break;
		  gout.write(buffer, 0, bytesRead2);
		}
		gout.close();
	      }
	      catch (IOException ioe){
		System.out.println("Error...URLHandler, getMatrixData");
	      }
	    }

	    else{
		FileWriter outMtrx = new FileWriter(data);
		int c;
		while (( c = in.read() ) != -1 )
		    outMtrx.write(c);
		in.close();
		outMtrx.close();
	    }
            
	}

	// If the URL is invalid, display an error (though it should be valid
        // if it reaches this point)
	catch(FileNotFoundException fnfe){
	    System.out.println("Filenotfound in URLExtractor.java");
	    return false;
	}

	// If there is no data on the URL's page, return false
	if (data == null || data == ""){
	    System.out.println("New false checker");
	    return false;
	}
	else
	    return true;
	
    } // End of method getMatrixData
    
} // end of public class URLHandler
    
