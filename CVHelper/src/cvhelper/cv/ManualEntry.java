package cvhelper.cv;

import bibtex.dom.*;

/**
 * ManualEntry
 * Extension of Entry class, implements special handling rules of ManualEntry BibTeX entry
 *
 */
public class ManualEntry extends Entry {

	/**
	 * ManualEntry
	 * Basic constructor, creates a new Entry using the underlying BibTeXEntry element.
	 *
	 */
	ManualEntry(BibtexEntry entry);

	/**
	 * ManualEntry
	 * Basic constructor, a kind of "copy" constructor between BibTeX types, converts the underlying element to this type
	 *
	 */
	public class ManualEntry(Entry entry);

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
	 * setOrganization
	 * Sets the value of the Organization BibTeX field in the BibTeX Entry
	 *
	 */
	public void setOrganization(String Organization);

	/**
	 * getOrganization
	 * Gets the current value of the Organization BibTeX field in the BibTeX Entry
	 *
	 */
	public String getOrganization();

	/**
	 * hasOrganization
	 * Checks whether the Organization BibTeX field in the BibTeX Entry has been initialized
	 *
	 */
	public boolean hasOrganization();

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
	 * setEdition
	 * Sets the value of the Edition BibTeX field in the BibTeX Entry
	 *
	 */
	public void setEdition(int Edition);

	/**
	 * getEdition
	 * Gets the current value of the Edition BibTeX field in the BibTeX Entry
	 *
	 */
	public int getEdition();

	/**
	 * hasEdition
	 * Checks whether the Edition BibTeX field in the BibTeX Entry has been initialized
	 *
	 */
	public boolean hasEdition();

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
