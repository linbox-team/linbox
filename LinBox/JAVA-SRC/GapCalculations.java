/**
 * @author Brian Gold
 * @version 1.1
 *
 * Date: 11/14/01
 * Last Updated: 1/7/02
 * File: GapCalculations.java
 *
 * This file contains some miscellaneous methods used by GapServlet and
 * SmithFormServlet in building the user-selected GAP commands
 */


public class GapCalculations{
    /**
     * @param none
     * @return StringBuffer
     *
     * This method constructs the correct GAP command based on the amount
     * of comments that the user wants and returns a StringBuffer with the
     * command
     */
    public static StringBuffer getCommentaryCommand(String comments){
	StringBuffer sbComments = new StringBuffer("");
    
	// If the user does not want any commands
	if (comments.equals("none") ){
	    sbComments.append("HomologyInfo(0);");
	}

	// If the user wants some commands
	else if (comments.equals("some") ){
	    sbComments.append("HomologyInfo(15);");
	}

	// If the user wants verbose comments
	else if (comments.equals("verbose") ){
	    sbComments.append("HomologyInfo(1000);");
	}

	// Just in case...should never be here
	else
	    System.out.println("In getCommentaryCommand()..shouldn't be here");
	
	// return the command as a StringBuffer
	return sbComments;
	
    }  // end of method getCommentaryCommand

    /**
     * @param none
     * @return StringBuffer
     *
     * This method constructs the correct GAP command based on the level of
     * timing that the user wants and returns a StringBuffer with the
     * command
     */
    public static StringBuffer getTimingCommand(String timing){
	StringBuffer sbTiming = new StringBuffer("");
	
	// if the user wants no timings
	if (timing.equals("none") ){
	    sbTiming.append("SetInfoLevel(InfoTiming, 0);");
	}

	// if the user wants high level timings
	else if (timing.equals("high level") ){
	    sbTiming.append("SetInfoLevel(InfoTiming, 10);");
	}

	// if the user wants detail timings
	else if (timing.equals("detail") ){
	    sbTiming.append("SetInfoLevel(InfoTiming, 100);");
	}

	// Just in case ... Should never be here
	else
	    System.out.println("In getTimingCommand()...shouldn't be here");
	
	// return the command as a StringBuffer
	return sbTiming;
	
    }  // end of method getTimingCommand


    /**
     * @param String matrixFilename
     * @param boolean boolMod
     * @param int intMod
     * @return StringBuffer
     *
     * This method determines which calculation the user wanted done
     * from the SmithForm choices and constructs the actual GAP
     * commands.  If the boolMod is false, the user has not entered
     * a modulus, and thus intMod is not even read.  Otherwise,
     * the modulus is used in construction of the command.  It returns
     * a StringBuffer containing the exact command to be run.  The string
     * matrixFilename is a string containing the name of the local file
     * the matrix is saved as, and is needed to construct the commands
     */
    public static StringBuffer getSmithCommand(String matrixFilename,
					 boolean boolMod, int intMod,
                                         String smithChoice, String algorithm){
	
	StringBuffer sbSmith = new StringBuffer("");
	
	// if the user wants the Smith Normal form (compact output form)
	if (smithChoice.equals("compact") ){
	    sbSmith.append("SMInvariantFactors(\"").append(matrixFilename);
	    sbSmith.append("\"");
	    if (boolMod == false){
		sbSmith.append(":HomologyAlgorithm:=ValenceElimAlgorithm);");
	    }
	    else{
		sbSmith.append(", ").append(intMod).append(":HomologyAlgorithm:=ValenceElimAlgorithm);");
	    }
	}

	// If the user wants the Smith Normal Form (long output form)
	else if (smithChoice.equals("long") ){
	    sbSmith.append("SMSmithForm(\"").append(matrixFilename);
	    sbSmith.append("\"");
	    if (boolMod == false){
		sbSmith.append(":HomologyAlgorithm:=ValenceElimAlgorithm);");
	    }
	    else{
		sbSmith.append(", ").append(intMod).append(":HomologyAlgorithm:=ValenceElimAlgorithm);");
	    }
	}

	// If the user wants the rank instead
	else if (smithChoice.equals("rank") ){
	    if (boolMod == false){
		sbSmith.append("SMIntegerRank(\"").append(matrixFilename);
		sbSmith.append("\"");
                if (algorithm.equals("elim") )
                    sbSmith.append(":ValuenceElimAlgorithm);");
                else
		    sbSmith.append(":HomologyAlgorithm:=ValenceBBAlgorithm);");
	    }
	    else{
		sbSmith.append("SMPrimePowerRank(\"").append(matrixFilename);
		sbSmith.append("\"");
                if (algorithm.equals("elim") )
                    sbSmith.append(", ").append(intMod).append(":HomologyAlgorithm:=ValenceElimAlgorithm);");
                else
		    sbSmith.append(", ").append(intMod).append(":HomologyAlgorithm:=ValenceBBAlgorithm);");
	    }
	}
	
	// Return the command as a StringBuffer
	return sbSmith;
	
    }  // end of method getSmithCommand



} // end of class GapCalculations
