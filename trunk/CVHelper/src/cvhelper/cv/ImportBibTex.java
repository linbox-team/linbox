package cvhelper.cv;


import java.io.*;
import java.util.*;
import bibtex.parser.BibtexParser;

public class ImportBibTex {

    public static List parseBibTeX(String bibtexFilename, CV cv) throws IOException {
        // 1.  Try to open the BibTeX filename
        // 2.  Use the BibtexParser to parse the BibTex
        // 3.  Use the CV's "makeNewEntry" Entry factory method to get
        // a new entry for each element parsed out of the imported bibtex file
        // 4.  Also use the CV's "getEntry" method to check to make sure there
        // aren't any dupes in the list
        // 5.  Return the list that first has all the entry's read out of the
        // user supplied BibTeX file.  The last two elements of this list will then be
        // Integers, where the first is the number of entries that had Parse errors
        // and were therefore not read, and the second was the list of entries that
        // were duplicates and were therefore not read.  This info is important
        // for the outputting of the Import BibTeX operation.
        // 6.  For any unrecoverable error that results in no BibTeX being read,
        // an IOException will be thrown.  For an unrecoverable error that prevents
        // further BibTeX from being read, the number indicating # of parse errored
        // entries will be returned negative (the negative indicating that further
        // content of the file could not be read, and was therefore abandoaned).

        return null;
    }

}
