package cvhelper;

import cvhelper.cv.*;
import cvhelper.ui.*;



/**
 *
 * <p>Title: CVHelperMain</p>
 *
 * <p>Description: The Main class that is executed at the start of the system (ie, to run
 *    this program from the command line, you would type java cvhelper.CVHelperMain. This
 *    Class figures out the user's home directory, finds the appplication data structure,
 *    creates the CV class, and initializes/launches any and all GUI threads</p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: </p>
 *
 * @author Rich Seagraves
 * @version 1.0
 */

public class CVHelperMain {
    public static void main(String[] args) {

        // this class does the following
        // 1.  Figures out the user's home directory and from that finds the BibTeX files and
        // config files.  If they don't exist, create them.
        // 2.  Create the CV class using Readers of these files
        // 3.  Initialize the different threads for the UI classes, and launch
        // the main UI.

    }
}
