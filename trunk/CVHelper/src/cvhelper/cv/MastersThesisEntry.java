package cvhelper.cv;

import bibtex.dom.*;

/**
 * MastersThesisEntry
 * Extension of Entry class, implements special handling rules of MastersThesisEntry BibTeX entry
 *
 */
public class MastersThesisEntry extends Entry {

	/**
	 * MastersThesisEntry
	 * Basic constructor, creates a new Entry using the underlying BibTeXEntry element.
	 *
	 */
	MastersThesisEntry(BibtexEntry entry);

	/**
	 * MastersThesisEntry
	 * Basic constructor, a kind of "copy" constructor between BibTeX types, converts the underlying element to this type
	 *
	 */
	public class MastersThesisEntry(Entry entry);

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
	 * setSchool
	 * Sets the value of the School BibTeX field in the BibTeX Entry
	 * Note that the field School is a required field of this BibTeX Type
	 *
	 */
	public void setSchool(String School);

	/**
	 * getSchool
	 * Gets the current value of the School BibTeX field in the BibTeX Entry
	 * Note the School field is a required field of this type of BibTeX entry
	 *
	 */
	public String getSchool();

	/**
	 * hasSchool
	 * Checks whether the School BibTeX field in the BibTeX Entry has been initialized
	 *
	 */
	public boolean hasSchool();

	/**
	 * setYear
	 * Sets the value of the Year BibTeX field in the BibTeX Entry
	 * Note that the field Year is a required field of this BibTeX Type
	 *
	 */
	public void setYear(int Year);

	/**
	 * getYear
	 * Gets the current value of the Year BibTeX field in the BibTeX Entry
	 * Note the Year field is a required field of this type of BibTeX entry
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
		return super.isValid() && hasTitle() && hasAuthor() && hasSchool() && hasYear();
	}

        // This method is a simple extension to Entry's clone method
        public Object clone() { return null;}

}
