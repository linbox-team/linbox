package cvhelper.cv;

import dom.*;
import java.util.Hashtable;
import java.util.Vector;


/**
 * <p>Title: </p>
 * <p>Description: </p>
 * This class manages all CV entries in the system.
* It has a collection of CV entries that have either been imported or manually input.
 * <p>Copyright: Copyright (c) 2004</p>
 * <p>Company: </p>
* Arrow Technologies
 * @author Ning Kang
 * @version 1.0
 */

public class CV {

  private BibtexFile bibfile;
  private Hashtable entryTable;   // mapping entry id to Entry object
  private Vector categories;


  protected CV() {
    /*
    bibfile = new BibtexFile() ;
    */
  }

  /**
   * CV
   * is a copy constructor
   * @param CV
   */
  public CV(CV newcv) {
    /*
    copy fields of newcv to field of cv.
    */
  }

  /**
   * CV(Reader bibReader, Reader configReader)
   * a constructor which contructs CV from BibTex file and configuration file.
   * @param bibReader
   * is an Reader object that reads BibTex file, from which bibliographical information
   * are putted into CV.
   *
   * @param configReader
   * is a Reader object that reads configuration file, from which categories and names
   * are putted into CV.
   */

  public CV(Reader bibReader, Reader configReader) {
    /*
    bibfile = new BibtexFile ;
    BibtexParser parser = new BibtexParser ;
    parser grabs information from bibReader into bibfile ;
    parser grabs information from configReader into bibfile;
    */
  }

  /**
   * createEntry
   * Creates a new entry object which is not put into CV yet.
   *
   * @param entryType
   * is the type of entry would be created
   *
   @ param entryId
   * is the bibtex id for the entry.
   *
   * @return:
   * return a new Entry with entryId.

   */

  public Entry createEntry(String entryType, String entryId){
    /*
    if ( entryId in the entryTable) return null ;
    else return new Entry with the entryId and entryType;
    */
  }

  /**
   * addEntry
   * adds an new entry into CV.
   * @param entry
   * An Entry object to be added into CV.
   *
   * @return:
   * return true if adding an entry is successful, else return false.
   */

  public boolean addEntry(Entry entry) {
    /*
    if (the id of the entry in entryTable ) return false ;
    else {
      add entry to bibfile;
      put entry to entryTable ;
      return true ;
    }
    */
  }

  /**
   * removeEntry
   * deletes an Entry object from CV.
   * @param entry
   * An Entry object need to be removed from CV.
   */
  public void removeEntry(Entry entry) {
    /*
    remove the entry from bibfile and entryTable;
    */
  }


  /**
   * getEntry
   * gets an entry by entryId.
   * @param entryId
   * bibtex id for a Entry object.
   *
   * @return:
   * return an Entry object if it is in the CV.
   */

  public Entry getEntry(String entryId) {
    /*
    if( entryTable contains the entryId key)
      return Entry with entryId ;
    else return null ;
    }
    */
  }

  /**
   * getEntries
   *
   * @return:
   * returns an array of entries in CV.
   */
  public Entry[] getEntries() {
    /*
    return entries in bibtexfile;
    */
  }



  /**
   * changeEntryId
   * changes the bibtex id of an entry.
   *
   * @ param newId
   * new bibtex id for the entry.
   *
   * @ param entry
   * the entry need to change id.

   * @return:
   * return a new Entry with entryId.

   */

  public boolean changeEntryId(String newId, Entry entry){
    /*
    if ( newId in the entryTable) return false ;
    else {
      change bibtex id of the entry to newId
      return true;
    }
    */
  }

  /**
   * getCategories
   * returns an array of categories stored in field categories.
   *
   * @return:
   * a list of categories
   */
  public String[] getCategories() {
    /*
    return an array of categories ;
    */
  }

  /**
   * setCategories
   * set a list of categories to the field of categories
   *
   *@param categories
   * a list of categories
   */

  private void setCategories(String[] categories) {
    /*
    set a list of categories to field of categories ;
    */
  }

  /**
   * addCategory
   * adds an category to the field of categories.
   *
   * @param category
   * the name of category
   *
   * @return:
   * return true if adding an category is successful, else return false.
   */
  public boolean addCategory(String category) {
    /*
    if ( category doesn't exist in categories) ) {
      add value category to field categories;
      return true ;
    }
    return false ;
    */
  }

  /**
   * removeCategory
   * removes an category.
   *
   * @param category
   * the name of category
   */
  public void removeCategory(String category) {
    /*
    if (category in the categories) {
      remove all the entries belong to that category ;
      remove the category from field categories ;
    }
    */
  }


  /**
   * changeCategory
   * changes old category to new category.
   *
   * @param oldone
   * the name of ole category
   *
   * @param newone
   * the name of new category
   * @return:
   * return true if changing an category is successful, else return false.
   */
  public boolean changeCategory(String oldone, String newone) {
    /*
    if ( category exists in categories) ) {
      change the old category of all the entries to new category;
      update the category in field categories ;
      return true ;
    }
    return false ;
    */
  }


  /**
   * getEntriesByCategories
   * gets a hashtable of categories to the list of entries. The category is the key.
   * List of entries is the mapping value.
   *
   * @return:
   * return true if changing an category is successful, else return false.
   */
  public Hashtable getEntriesByCategories() {
    /*
    create a hashtable ;
    for ( category in the categories )
    {
       find all the entries in the CV belongs to the category ;
       insert the pair of category and the list of entries into the hashtable ;
    }
    return the hashtable ;
    */
  }

  /**
   * writeCV
   * writes content of CV using writer.
   *
   */
  public void writeCV(PrintWriter writer) {
    /*
     writes out all entries of CV using writer;
    */
  }

  public void writeConfig(PrintWriter writer) {
   /*
     writes out categories and names using writer;
   */
  }

  /**
   * createBibtexPerson
   * creates a BibtexPerson object. This function has package access limitation.
   * Only the classes in the package can call this function.
   *
   */
  BibtexPerson createBibtexPerson(String first,
        String preLast,
        String last,
        String lineage,
        boolean isOthers
) {
/*
    return a new BibtexPerson with given parameters.
 */
  }

  /**
   * createPersonList
   * creates a createPersonList object. This function has package access limitation.
   * Only the classes in the package can call this function.
   *
   */

  public BibtexPersonList createPersonList(){
      return new BibtexPersonList;
  }


  /**
   * createString
   * creates a BibtexString object. This function has package access limitation.
   * Only the classes in the package can call this function.
   * Content does not include the quotes or curly braces around the string!
   *
   * @param content
   * string need to create a BibtexString
   */


        public BibtexString createString(String content){
                return new BibtexString;
        }

  /**
   * Templates and Filters functions will be handled in next phase.
   *
   */

}
