package cvhelper.cv;

import bibtex.dom.*;

/**
 * UnPublishedEntry
 * Extension of Entry class, implements special handling rules of UnPublishedEntry BibTeX entry
 *
 */
public class UnPublishedEntry extends Entry {

	/**
	 * UnPublishedEntry
	 * Basic constructor, creates a new Entry using the underlying BibTeXEntry element.
	 *
	 */
	UnPublishedEntry(BibtexEntry entry);

	/**
	 * UnPublishedEntry
	 * Basic constructor, a kind of "copy" constructor between BibTeX types, converts the underlying element to this type
	 *
	 */
	public class UnPublishedEntry(Entry entry);

	/**
	 * setTitle
	 * Sets the value of the Title BibTeX field in the BibTeX Entry
	 * Note that the field Title is a required field of this BibTeX Type
	 *
	 */
	public void setTitle(String Title);

	/**
	 * getTitle
	 * Gets the current value of the Title BibTeX field in the BibTeX Entry
	 * Note the Title field is a required field of this type of BibTeX entry
	 *
	 */
	public String getTitle();

	/**
	 * hasTitle
	 * Checks whether the Title BibTeX field in the BibTeX Entry has been initialized
	 *
	 */
	public boolean hasTitle();

	/**
	 * setAuthor
	 * Sets the value of the Author BibTeX field in the BibTeX Entry
	 * Note that the field Author is a required field of this BibTeX Type
	 *
	 */
	public void setAuthor(List L);

	/**
	 * getAuthor
	 * Gets the current value of the Author BibTeX field in the BibTeX Entry
	 * Note the Author field is a required field of this type of BibTeX entry
	 *
	 */
	public List getAuthor();

	/**
	 * hasAuthor
	 * Checks whether the Author BibTeX field in the BibTeX Entry has been initialized
	 *
	 */
	public boolean hasAuthor();

	/**
	 * setNote
	 * Sets the value of the Note BibTeX field in the BibTeX Entry
	 * Note that the field Note is a required field of this BibTeX Type
	 *
	 */
	public void setNote(String Note);

	/**
	 * getNote
	 * Gets the current value of the Note BibTeX field in the BibTeX Entry
	 * Note the Note field is a required field of this type of BibTeX entry
	 *
	 */
	public String getNote();

	/**
	 * hasNote
	 * Checks whether the Note BibTeX field in the BibTeX Entry has been initialized
	 *
	 */
	public boolean hasNote();

	/**
	 * setMonth
	 * Sets the value of the Month BibTeX field in the BibTeX Entry
	 *
	 */
	public void setMonth(String Month);

	/**
	 * getMonth
	 * Gets the current value of the Month BibTeX field in the BibTeX Entry
	 *
	 */
	public String getMonth();

	/**
	 * hasMonth
	 * Checks whether the Month BibTeX field in the BibTeX Entry has been initialized
	 *
	 */
	public boolean hasMonth();

	/**
	 * setYear
	 * Sets the value of the Year BibTeX field in the BibTeX Entry
	 *
	 */
	public void setYear(int Year);

	/**
	 * getYear
	 * Gets the current value of the Year BibTeX field in the BibTeX Entry
	 *
	 */
	public int getYear();

	/**
	 * hasYear
	 * Checks whether the Year BibTeX field in the BibTeX Entry has been initialized
	 *
	 */
	public boolean hasYear();

	/**
	 * isValid
	 * Checks whether all the requried fields of this entry have been defined using the has{Field} methods.
	 *
	 */
	public boolean isValid() {
		return super.isValid() && hasTitle() && hasAuthor() && hasNote();
	}

        // This method is a simple extension to Entry's clone method
        public Object clone() { return null;}

}
