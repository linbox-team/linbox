import java.net.*;
import java.io.*;


/**
 * @author Brian Gold
 * @version 1.0
 * 
 * Date: 7/25/01
 * Last Updated: 11/14/01
 * File: SendMail.java
 * 
 * This file contains only static methods which are used to attempt to
 * send e-mail notifying the user that the GAP/Homology server has finished
 * the calculations and that the results can be viewed now, along with the
 * random number used to retrieve the results.
 */
public class SendMail{

    /**
     * @param String to
     * @param String from
     * @param String subj
     * @param String data
     * @return none
     *
     * This method sends the mail to the user whose email address is contained
     * in the string 'to', informing the user that the calculations have
     * completed.  This method is called by GapServlet.
     */
    public static synchronized void sendMail(String to, String from,
					     String subj, String data)
	throws IOException, ProtocolException{

	try {
	    // Try to establish a connection to port 25 on mail.eecis.udel.edu
	    // which is the port used for sending e-mail
	    //System.out.println("Establishing connection..");
	    Socket sock=new Socket("mail.eecis.udel.edu", 25);

	    // Set up both an input and an output stream for the server
	    // to send data out and to read any replies in
	    DataOutputStream out=new DataOutputStream(sock.getOutputStream());
	    DataInputStream in=new DataInputStream(sock.getInputStream());
	    //System.out.println("Connection OK, sending data");
	    
	    if(readReply(in)) {
		sock.close();
		return;
	    }

	    // Send a greeting to the mail server and read the server's reply
	    out.writeBytes("HELO mail.eecis.udel.edu\r\n");
	    if(readReply(in)) {
		sock.close();
		return;
	    }

	    // If there are no errors, send the email address of the
	    // sender of the email and read the server's reply
	    out.writeBytes("MAIL FROM:<"+from+">\r\n");
	    if(readReply(in)){
		sock.close();
		return;
	    }

	    // If there are no errors, send the email address of the
	    // recipient of the email and read the server's reply
	    out.writeBytes("RCPT TO:<"+to+">\r\n");
	    if(readReply(in)) {
		sock.close();
		return;
	    }

	    // If there are no errors, send the data of the email
	    // and read the server's reply
	    out.writeBytes("DATA\r\n");
	    if(readReply(in)) {
		sock.close();
		return;
	    }

	    // If there is a subject field, send it, and read
	    // the server's reply
	    if(subj!=null)
		out.writeBytes("Subject: "+subj+"\n\n");
	    out.writeBytes(data+"\r\n.\r\nQUIT\r\n");
	    if(readReply(in)) {
		sock.close();
		return;
	    }

	    // Finish sending the mail and close the port
	    sock.close();

	} // end of try block

	// If there is a problem setting up the port, display an error
	// to the server (not to the user, though this might be amended later)
	catch(IOException e) {
	    System.out.println("Connection problem or connection closed");
	}

    }  // end of method sendMail


    /**
     * @param DataInputStream in
     * @return boolean
     * 
     * This method is called by sendMail() and determines if there are
     * any errors in the sending process by reading the reply from the
     * server.  If there are errors, this method
     * returns TRUE.  If there are NO errors, this method returns FALSE
     */
    public static boolean readReply(DataInputStream in) {

	String line= "";

	// Read the reply from the server
	try {
	    line=in.readLine();
	}

	// There is an error reading the reply from the server; mail
	// will not be send; return true
	catch(IOException e) {
	    System.out.println("Reading reply problem");
	    return true;
	}

	// If the server reply starts with a 4 or a 5, there is an error
	// so return true
	if(line.startsWith("4") || line.startsWith("5")) {
	    System.out.println("Err: "+line);
	    return true;
	}

	// Else, there is no error, so return false so sendMail() can
	// continue trying to send the message to the user
	else 
	    return false;

    }  // end of method readReply


}  // end of public class SendMail
