package mypackage;

import java.util.*;
import java.util.zip.*;
import java.net.*;
import java.io.*;
import java.text.*;
import java.sql.*;
import javax.servlet.*;
import javax.servlet.http.*;



/**
 * @author Eric Schrag
 * @version 1.2
 * 
 * Date: 6/13/02
 * Last Updated: 1/07/03
 * File: serverSearch.java
 * 
 * This class connects to a remote database, which will be in charge of a 
 * list of servers available and whether or not they are in use, and uses
 * the following methods:
 * getServer: retrieves an available server name
 * changeStatusUsed: changes an IP address' status to in use in the database
 * checkServer: boolean function, returns true if a server is available
 */


public class serverSearch{

    public Connection conn;
    private Statement stmt;
    private ResultSet rs;
    private String addr;


    public serverSearch(){

	System.out.println("Checking for server");

    this.stmt = null;
    this.rs = null;
    this.addr = null;
    this.conn = null;


    try{

    Class.forName("com.mysql.jdbc.Driver");
    this.conn = DriverManager.getConnection("jdbc:mysql://linbox.pc.cis.udel.edu:3306/linbox?user=root&autoReconnect=true");

    }

	catch(ClassNotFoundException E) {
	        System.out.println("Error in ClassCreation of driver");
	}

	catch(SQLException E) {
	    System.out.println("Error in SQL");
	    System.out.println(E.getMessage());
	    System.out.println(E.getSQLState());
	    E.printStackTrace();
	}

	catch(Exception e) {
	                System.out.println("General Error");
	    }

    }



    public String getServer(){


	System.out.println("Attempting to get server");

	System.out.println(this.conn);

	this.addr = null;

	try{


	    //	    if (this.conn != null){

	       this.stmt = this.conn.createStatement();

	       this.rs = stmt.executeQuery("SELECT ip FROM servers WHERE in_use = 0 AND down = 0");

	       this.rs.next();

	       this.addr = this.rs.getString("ip");

	       return (this.addr);
	       //	    }



	}



	catch(SQLException E) {
	    System.out.println("Debug A");
	}

	catch(Exception e) {
	    System.out.println("Debug B");
	    System.out.println(e);
	    }

	finally{

	    return(this.addr);

	}



    }


    public void changeStatusUsed(){

	try{



	       this.stmt = this.conn.createStatement();

	       stmt.executeQuery("UPDATE servers SET in_use=1 WHERE ip = '" + InetAddress.getLocalHost().getHostName() + "'");



	}



	catch(SQLException E) {
	    System.out.println("Debug A");
	}

	catch(Exception e) {
	    System.out.println("Debug B");
	    }

    }




    public Boolean checkServer(){


	System.out.println("Attempting to get server");

	System.out.println(this.conn);

	this.addr = null;

	try{


	    //	    if (this.conn != null){

	       this.stmt = this.conn.createStatement();

	       this.rs = stmt.executeQuery("SELECT ip FROM servers WHERE in_use = 0 AND down = 0");

	       this.rs.next();

	       this.addr = this.rs.getString("ip");

	       if(this.addr == null)
		   return(Boolean.FALSE);
	       else
		   return(Boolean.TRUE);
	       //	    }



	}



	catch(SQLException E) {
	    System.out.println("Debug A");
	}

	catch(Exception e) {
	    System.out.println("Debug B");
	    System.out.println(e);
	    }

	finally{

	       if(this.addr == null)
		   return(Boolean.FALSE);
	       else
		   return(Boolean.TRUE);

	}



    }



}
