import java.util.*;
import java.util.zip.*;
import java.net.*;
import java.io.*;
import java.text.*;
import java.lang.Thread;
import java.sql.*;
import javax.servlet.*;
import javax.servlet.http.*;


/*
 * class: serverDaemon
 * version: 0.9
 * Programmer: Eric Schrag
 *
 * Purpose: This is java that can run quietly in the background
 *          in a loop, accessing a database to get a server list
 *          and checking each of those servers to see if they are up.
 *          if they are not, the class will set the servers "down" flag
 *          to 1, in order to keep LinBox's server from thinking the
 *          server is available for processing. If they are available,
 *          no change is made in the database other than to mark "down"
 *          as 0.
 **/




public class serverDaemon{

    public static void main(String[] stuff){
	
	Process proc;
	Connection conn;
        Statement stmt;
	ResultSet rs;
	String addr;
	Boolean connflag;
	

	//	connflag = Boolean.TRUE;
	
	while(true){
	    
	    try{
		Thread.sleep(10000);
	    }
	    catch(Exception e) {
		System.out.println("General Error");
	    }
	    

	    
	    
	    //	    System.out.println("Checking for server");
	    
	    stmt = null;
	    rs = null;
	    addr = null;
	    conn = null;
	    
	    
	    try{
		
		Class.forName("com.mysql.jdbc.Driver");

	System.out.println("Hey, i do something, at least");
		
		connflag = Boolean.TRUE;

		conn = DriverManager.getConnection("jdbc:mysql://linbox.pc.cis.udel.edu/linbox?user=root");
	    }
	    
	    catch(ClassNotFoundException E) {
		System.out.println("Error in ClassCreation of driver");
		System.out.println(E.getMessage());
	    }
	    
	    catch(SQLException E) {
		connflag = Boolean.FALSE;
		System.out.println(E.getMessage());
	    }
	    
	    catch(Exception e) {
		System.out.println("General Error");
	    }

	    
	    
	    
	    //    System.out.println("Attempting to get server");
	    
	    // System.out.println(conn);
	    
	    addr = null;
	    
	    try{
		
		
		//		if( !connflag.booleanValue() ){
		//    System.out.println("MySQL server is down- restarting");
		//    proc = Runtime.getRuntime().exec("/usr/bin/safe_mysqld &");
		//    System.out.println("ooga booga");
		//}		
		
		stmt = conn.createStatement();
		
		rs = stmt.executeQuery("SELECT ip FROM servers");
		
		URL url;
		URLConnection siteCon;
		
	       	rs.first();
		while( !rs.isAfterLast() ){
		    
		    try{
			url = new URL( "http://" + rs.getString("ip") + ":8080");
		    }
		    catch(MalformedURLException murle){
			//	System.out.println(rs.getString("ip") + "was not a valid server");
			url = null;
		    }
		    
		    if(url != null){
			siteCon = url.openConnection();
		  
		    
			if(siteCon.getContentLength() <= 0 ){
			    stmt = conn.createStatement();

			    System.out.println("Server not found: " + rs.getString("ip") );
			    
			    if(rs.getString("ip").equals("linbox.pc.cis.udel.edu")){
			    System.out.println("restarting linbox");
				proc = Runtime.getRuntime().exec("startup.sh");
				proc.waitFor();
			    }
			    else{
			    stmt.executeQuery("UPDATE servers SET down=1 WHERE ip = '" + rs.getString("ip") + "'");
			    System.out.println(rs.getString("ip") + " is down");
			    }
			}
			else{
			    System.out.println(rs.getString("ip") + " is up");

			    stmt = conn.createStatement();
			    stmt.executeQuery("UPDATE servers SET down=0 WHERE ip = '" + rs.getString("ip") + "'");

			}


		    }

		    rs.next();
		}
		//	    }
		
		
		
	    }
	    
	    

	    catch(SQLException E) {
		//		System.out.println("Debug A");
		//System.out.println(E);
	    }
	    
	    catch(Exception e) {
		//System.out.println("Debug B");
		//System.out.println(e);
	    }
	    
	    
	}
	
    }
    
}
