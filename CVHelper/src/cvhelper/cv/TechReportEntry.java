package cvhelper.cv;

import bibtex.dom.*;

/**
 * TechReportEntry
 * Extension of Entry class, implements special handling rules of TechReportEntry BibTeX entry
 *
 */
public class TechReportEntry extends Entry {

	/**
	 * TechReportEntry
	 * Basic constructor, creates a new Entry using the underlying BibTeXEntry element.
	 *
	 */
	TechReportEntry(BibtexEntry entry);

	/**
	 * TechReportEntry
	 * Basic constructor, a kind of "copy" constructor between BibTeX types, converts the underlying element to this type
	 *
	 */
	public class TechReportEntry(Entry entry);

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
	 * setInstitution
	 * Sets the value of the Institution BibTeX field in the BibTeX Entry
	 * Note that the field Institution is a required field of this BibTeX Type
	 *
	 */
	public void setInstitution(String Institution);

	/**
	 * getInstitution
	 * Gets the current value of the Institution BibTeX field in the BibTeX Entry
	 * Note the Institution field is a required field of this type of BibTeX entry
	 *
	 */
	public String getInstitution();

	/**
	 * hasInstitution
	 * Checks whether the Institution BibTeX field in the BibTeX Entry has been initialized
	 *
	 */
	public boolean hasInstitution();

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
	 * setType
	 * Sets the value of the Type BibTeX field in the BibTeX Entry
	 *
	 */
	public void setType(String Type);

	/**
	 * getType
	 * Gets the current value of the Type BibTeX field in the BibTeX Entry
	 *
	 */
	public String getType();

	/**
	 * hasType
	 * Checks whether the Type BibTeX field in the BibTeX Entry has been initialized
	 *
	 */
	public boolean hasType();

	/**
	 * setNumber
	 * Sets the value of the Number BibTeX field in the BibTeX Entry
	 *
	 */
	public void setNumber(int Number);

	/**
	 * getNumber
	 * Gets the current value of the Number BibTeX field in the BibTeX Entry
	 *
	 */
	public int getNumber();

	/**
	 * hasNumber
	 * Checks whether the Number BibTeX field in the BibTeX Entry has been initialized
	 *
	 */
	public boolean hasNumber();

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
		return super.isValid() && hasTitle() && hasYear() && hasInstitution() && hasAuthor();
	}

        // This method is a simple extension to Entry's clone method
        public Object clone() { return null;}

}
