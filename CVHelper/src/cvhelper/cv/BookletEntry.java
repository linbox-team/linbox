package cvhelper.cv;

import bibtex.dom.*;

/**
 * BookletEntry
 * Extension of Entry class, implements special handling rules of BookletEntry BibTeX entry
 *
 */
public class BookletEntry extends Entry {

	/**
	 * BookletEntry
	 * Basic constructor, creates a new Entry using the underlying BibTeXEntry element.
	 *
	 */
	BookletEntry(BibtexEntry entry);

	/**
	 * BookletEntry
	 * Basic constructor, a kind of "copy" constructor between BibTeX types, converts the underlying element to this type
	 *
	 */
	public class BookletEntry(Entry entry);

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
	 *
	 */
	public void setAuthor(List L);

	/**
	 * getAuthor
	 * Gets the current value of the Author BibTeX field in the BibTeX Entry
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
	 * setHowPublished
	 * Sets the value of the HowPublished BibTeX field in the BibTeX Entry
	 *
	 */
	public void setHowPublished(String HowPublished);

	/**
	 * getHowPublished
	 * Gets the current value of the HowPublished BibTeX field in the BibTeX Entry
	 *
	 */
	public String getHowPublished();

	/**
	 * hasHowPublished
	 * Checks whether the HowPublished BibTeX field in the BibTeX Entry has been initialized
	 *
	 */
	public boolean hasHowPublished();

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
	 * setAddress
	 * Sets the value of the Address BibTeX field in the BibTeX Entry
	 *
	 */
	public void setAddress(String Address);

	/**
	 * getAddress
	 * Gets the current value of the Address BibTeX field in the BibTeX Entry
	 *
	 */
	public String getAddress();

	/**
	 * hasAddress
	 * Checks whether the Address BibTeX field in the BibTeX Entry has been initialized
	 *
	 */
	public boolean hasAddress();

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
	 * setNote
	 * Sets the value of the Note BibTeX field in the BibTeX Entry
	 *
	 */
	public void setNote(String Note);

	/**
	 * getNote
	 * Gets the current value of the Note BibTeX field in the BibTeX Entry
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
	 * isValid
	 * Checks whether all the requried fields of this entry have been defined using the has{Field} methods.
	 *
	 */
	public boolean isValid() {
		return super.isValid() && hasTitle();
	}

        // This method is a simple extension to Entry's clone method
        public Object clone() { return null;}


}
