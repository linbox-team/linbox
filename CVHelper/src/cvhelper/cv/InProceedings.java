package cvhelper.cv;

import bibtex.dom.*;

/**
 * InProceedings
 * Extension of Entry class, implements special handling rules of InProceedings BibTeX entry
 *
 */
public class InProceedingsEntry extends Entry {

	/**
	 * InProceedings
	 * Basic constructor, creates a new Entry using the underlying BibTeXEntry element.
	 *
	 */
	InProceedingsEntry(BibtexEntry entry);

	/**
	 * InProceedings
	 * Basic constructor, a kind of "copy" constructor between BibTeX types, converts the underlying element to this type
	 *
	 */
	public class InProceedingsEntry(Entry entry);

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
	 * setBookTitle
	 * Sets the value of the BookTitle BibTeX field in the BibTeX Entry
	 * Note that the field BookTitle is a required field of this BibTeX Type
	 *
	 */
	public void setBookTitle(String BookTitle);

	/**
	 * getBookTitle
	 * Gets the current value of the BookTitle BibTeX field in the BibTeX Entry
	 * Note the BookTitle field is a required field of this type of BibTeX entry
	 *
	 */
	public String getBookTitle();

	/**
	 * hasBookTitle
	 * Checks whether the BookTitle BibTeX field in the BibTeX Entry has been initialized
	 *
	 */
	public boolean hasBookTitle();

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
	 * setEditor
	 * Sets the value of the Editor BibTeX field in the BibTeX Entry
	 *
	 */
	public void setEditor(List L);

	/**
	 * getEditor
	 * Gets the current value of the Editor BibTeX field in the BibTeX Entry
	 *
	 */
	public List getEditor();

	/**
	 * hasEditor
	 * Checks whether the Editor BibTeX field in the BibTeX Entry has been initialized
	 *
	 */
	public boolean hasEditor();

	/**
	 * setVolumeOrNumber
	 * Sets the value of the VolumeOrNumber BibTeX field in the BibTeX Entry
	 *
	 */
	public void setVolumeOrNumber(int VolumeOrNumber);

	/**
	 * getVolumeOrNumber
	 * Gets the current value of the VolumeOrNumber BibTeX field in the BibTeX Entry
	 *
	 */
	public int getVolumeOrNumber();

	/**
	 * hasVolumeOrNumber
	 * Checks whether the VolumeOrNumber BibTeX field in the BibTeX Entry has been initialized
	 *
	 */
	public boolean hasVolumeOrNumber();

	/**
	 * setSeries
	 * Sets the value of the Series BibTeX field in the BibTeX Entry
	 *
	 */
	public void setSeries(String Series);

	/**
	 * getSeries
	 * Gets the current value of the Series BibTeX field in the BibTeX Entry
	 *
	 */
	public String getSeries();

	/**
	 * hasSeries
	 * Checks whether the Series BibTeX field in the BibTeX Entry has been initialized
	 *
	 */
	public boolean hasSeries();

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
	 * setPublisher
	 * Sets the value of the Publisher BibTeX field in the BibTeX Entry
	 *
	 */
	public void setPublisher(String Publisher);

	/**
	 * getPublisher
	 * Gets the current value of the Publisher BibTeX field in the BibTeX Entry
	 *
	 */
	public String getPublisher();

	/**
	 * hasPublisher
	 * Checks whether the Publisher BibTeX field in the BibTeX Entry has been initialized
	 *
	 */
	public boolean hasPublisher();

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
		return super.isValid() && hasTitle() && hasAuthor() && hasBookTitle() && hasYear();
	}

        // This method is a simple extension to Entry's clone method
        public Object clone() { return null;}

}
