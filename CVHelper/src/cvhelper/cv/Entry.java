package cvhelper.cv;

import bibtex.dom.BibtexEntry;


/**
 * <p>Title: Entry</p>
 *
 * <p>Description: The base class of the hierarchy of entries.  Provides
 *   methods that all entries should have.</p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: </p>
 *
 * @author Rich Seagraves
 * @version 1.0
 */
public class Entry {

    // This element is made package acessible so that
    // 1.  The subclasses of this class can access their own entry
    // 2.  The CV and ImportBibTex (especially ImportBibTex)
    // classes can access this field for special manipulations
    BibtexEntry _entry;

    // This is here for one reason:  to check to make sure when the user
    // tries to change the BibTeX ID, the entry first checks with the system
    // to ensure that this won't create multiple versions of the ID in the system
    private CV _cv;

    /**
     * Entry - Internal constructor, only meant to be used
     * by domain elements (CV, BibTeXimporter, etc)
     * Every Entry MUST have the BibtexEntry it is representing
     * This BibtexEntry is used to make the underlying changes to the
     * BibTeX
     * In addition, entries must have the CV from which it came, for
     * ID consistancy checking
     *
     * @param entry
     * @param cv
     */

  Entry(BibtexEntry entry, CV cv);


  /**
   * Entry - This constructor is meant to act as a copy
   * constructor and to be used by the UI code to convert
   * Entries from one type to another (ie, if you have
   * BookEntry e, you can automatically change it to an
   * ArticleEntry by saying
   * ArticleEntry ae = new ArticleEntry(e);
   * This will preserve all the necessary fields of the
   * particular entry and strip out the other elements of the
   * entry
   *
   * @param entry Entry
   */
  public Entry(Entry entry);

  //note that all of these calls are completed by interfacing through
  // the BibtexEntry (ie, either accessing it's underlying contents
  // or changing it's underlying contents

  // sets to one of the 14 types of BibTeX types.  Note that the
  // specialized type will attempt to correct and convert back to this type
  // (we hope this won't cause breakage, we'll see)
  public String getType();

  // expect one of the 14 types of BibTeX entries:
  // gives the type of the entry
  public void setType(String T);

  // returns true if the entry has a BibTex type
  // Mostly, this just means it is a defined BibTex file
  public boolean hasBibTexType() {
    return _entry != null;
  }

  // returns the category of this entry, which is stored in the
  // BibTeX representation.  Note - returns null if
  // hasCategory() == false
  public String getCategory();

  // expects some category of the categories defined in the user's
  // category list
  // Note that the category is written as an ignored field of each entry
  // (in the field CATEGORY = "User Defined Category")
  public void setCategory(String category);

  // checks to see whether a category has been defined for this
  // entry.  If not, return false
  public boolean hasCategory();

  // returns the bibTeXId, or citation key, of the element
  public String getID();

  // resets the citation key of the element to something else
  public boolean setID(String bibTeXID) {
      if(bibTeXID == getID())
          return true;
      else if(_cv.getEntry(bibTeXID) != null)
          return false;
      else {
          // change the ID first
          return true;
      }
  }

  // checks whether this thing has a BibTeXID
  public boolean hasID();

  /**
   * isValid
   * Checks to see whether the necessary fields of this particular
   * entry are valid. Notice that for the generic entry, all that is
   * needed is a category, a BibTeXId, and a BibTex entry type
   * Note that all subclasses of this class call this super class.
   *
   *
   * @return boolean
   */

  public boolean isValid() {
    return hasCategory() && hasID() && hasBibTexType();
  }


  // this method returns a copy of the entry, by creating a new entry that
  // contains all the new elements of the old entry
  // This method is the meat behind the clone methods of all the sub-classes
  // to Entry
  public Object clone() {return null;}


  /**
   * equals
   * Checks whether this entry is the same as another
   * Uses only the BibTeXID as the  distinguishing factor
   *
   * @param anObject Object
   * @return boolean
   */

  public boolean equals(Object anObject) {
      if(anObject instanceof Entry) {
          Entry otherEntry = (Entry) anObject;
          return getID() == otherEntry.getID();
      }
      return false;
  }

}
