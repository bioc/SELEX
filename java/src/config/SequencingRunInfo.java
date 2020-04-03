//
// This file was generated by the JavaTM Architecture for XML Binding(JAXB) Reference Implementation, v2.2.8-b130911.1802 
// See <a href="http://java.sun.com/xml/jaxb">http://java.sun.com/xml/jaxb</a> 
// Any modifications to this file will be lost upon recompilation of the source schema. 
// Generated on: 2020.04.02 at 08:24:17 PM PDT 
//


package config;

import java.util.ArrayList;
import java.util.List;
import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlType;


/**
 * <p>Java class for anonymous complex type.
 * 
 * <p>The following schema fragment specifies the expected content contained within this class.
 * 
 * <pre>
 * &lt;complexType>
 *   &lt;complexContent>
 *     &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
 *       &lt;sequence>
 *         &lt;element ref="{}DataFile"/>
 *         &lt;element ref="{}DataFileType" minOccurs="0"/>
 *         &lt;element ref="{}SequencingPlatform" minOccurs="0"/>
 *         &lt;element ref="{}ResearcherName" minOccurs="0"/>
 *         &lt;element ref="{}ResearcherEmail" minOccurs="0"/>
 *         &lt;element ref="{}SequencingFacilityName" minOccurs="0"/>
 *         &lt;element ref="{}SequencingFacilityEmail" minOccurs="0"/>
 *         &lt;element ref="{}Description" minOccurs="0"/>
 *         &lt;element ref="{}Notes" minOccurs="0"/>
 *         &lt;element ref="{}Sample" maxOccurs="unbounded"/>
 *       &lt;/sequence>
 *       &lt;attribute name="name" use="required" type="{http://www.w3.org/2001/XMLSchema}string" />
 *     &lt;/restriction>
 *   &lt;/complexContent>
 * &lt;/complexType>
 * </pre>
 * 
 * 
 */
@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "", propOrder = {
    "dataFile",
    "dataFileType",
    "sequencingPlatform",
    "researcherName",
    "researcherEmail",
    "sequencingFacilityName",
    "sequencingFacilityEmail",
    "description",
    "notes",
    "sample"
})
@XmlRootElement(name = "SequencingRunInfo")
public class SequencingRunInfo {

    @XmlElement(name = "DataFile", required = true)
    protected String dataFile;
    @XmlElement(name = "DataFileType", defaultValue = "FASTQ.TXT")
    protected String dataFileType;
    @XmlElement(name = "SequencingPlatform", defaultValue = "")
    protected String sequencingPlatform;
    @XmlElement(name = "ResearcherName", defaultValue = "")
    protected String researcherName;
    @XmlElement(name = "ResearcherEmail", defaultValue = "")
    protected String researcherEmail;
    @XmlElement(name = "SequencingFacilityName", defaultValue = "")
    protected String sequencingFacilityName;
    @XmlElement(name = "SequencingFacilityEmail", defaultValue = "")
    protected String sequencingFacilityEmail;
    @XmlElement(name = "Description", defaultValue = "")
    protected String description;
    @XmlElement(name = "Notes", defaultValue = "")
    protected String notes;
    @XmlElement(name = "Sample", required = true)
    protected List<Sample> sample;
    @XmlAttribute(name = "name", required = true)
    protected String name;

    /**
     * Gets the value of the dataFile property.
     * 
     * @return
     *     possible object is
     *     {@link String }
     *     
     */
    public String getDataFile() {
        return dataFile;
    }

    /**
     * Sets the value of the dataFile property.
     * 
     * @param value
     *     allowed object is
     *     {@link String }
     *     
     */
    public void setDataFile(String value) {
        this.dataFile = value;
    }

    /**
     * Gets the value of the dataFileType property.
     * 
     * @return
     *     possible object is
     *     {@link String }
     *     
     */
    public String getDataFileType() {
        return dataFileType;
    }

    /**
     * Sets the value of the dataFileType property.
     * 
     * @param value
     *     allowed object is
     *     {@link String }
     *     
     */
    public void setDataFileType(String value) {
        this.dataFileType = value;
    }

    /**
     * Gets the value of the sequencingPlatform property.
     * 
     * @return
     *     possible object is
     *     {@link String }
     *     
     */
    public String getSequencingPlatform() {
        return sequencingPlatform;
    }

    /**
     * Sets the value of the sequencingPlatform property.
     * 
     * @param value
     *     allowed object is
     *     {@link String }
     *     
     */
    public void setSequencingPlatform(String value) {
        this.sequencingPlatform = value;
    }

    /**
     * Gets the value of the researcherName property.
     * 
     * @return
     *     possible object is
     *     {@link String }
     *     
     */
    public String getResearcherName() {
        return researcherName;
    }

    /**
     * Sets the value of the researcherName property.
     * 
     * @param value
     *     allowed object is
     *     {@link String }
     *     
     */
    public void setResearcherName(String value) {
        this.researcherName = value;
    }

    /**
     * Gets the value of the researcherEmail property.
     * 
     * @return
     *     possible object is
     *     {@link String }
     *     
     */
    public String getResearcherEmail() {
        return researcherEmail;
    }

    /**
     * Sets the value of the researcherEmail property.
     * 
     * @param value
     *     allowed object is
     *     {@link String }
     *     
     */
    public void setResearcherEmail(String value) {
        this.researcherEmail = value;
    }

    /**
     * Gets the value of the sequencingFacilityName property.
     * 
     * @return
     *     possible object is
     *     {@link String }
     *     
     */
    public String getSequencingFacilityName() {
        return sequencingFacilityName;
    }

    /**
     * Sets the value of the sequencingFacilityName property.
     * 
     * @param value
     *     allowed object is
     *     {@link String }
     *     
     */
    public void setSequencingFacilityName(String value) {
        this.sequencingFacilityName = value;
    }

    /**
     * Gets the value of the sequencingFacilityEmail property.
     * 
     * @return
     *     possible object is
     *     {@link String }
     *     
     */
    public String getSequencingFacilityEmail() {
        return sequencingFacilityEmail;
    }

    /**
     * Sets the value of the sequencingFacilityEmail property.
     * 
     * @param value
     *     allowed object is
     *     {@link String }
     *     
     */
    public void setSequencingFacilityEmail(String value) {
        this.sequencingFacilityEmail = value;
    }

    /**
     * Gets the value of the description property.
     * 
     * @return
     *     possible object is
     *     {@link String }
     *     
     */
    public String getDescription() {
        return description;
    }

    /**
     * Sets the value of the description property.
     * 
     * @param value
     *     allowed object is
     *     {@link String }
     *     
     */
    public void setDescription(String value) {
        this.description = value;
    }

    /**
     * Gets the value of the notes property.
     * 
     * @return
     *     possible object is
     *     {@link String }
     *     
     */
    public String getNotes() {
        return notes;
    }

    /**
     * Sets the value of the notes property.
     * 
     * @param value
     *     allowed object is
     *     {@link String }
     *     
     */
    public void setNotes(String value) {
        this.notes = value;
    }

    /**
     * Gets the value of the sample property.
     * 
     * <p>
     * This accessor method returns a reference to the live list,
     * not a snapshot. Therefore any modification you make to the
     * returned list will be present inside the JAXB object.
     * This is why there is not a <CODE>set</CODE> method for the sample property.
     * 
     * <p>
     * For example, to add a new item, do as follows:
     * <pre>
     *    getSample().add(newItem);
     * </pre>
     * 
     * 
     * <p>
     * Objects of the following type(s) are allowed in the list
     * {@link Sample }
     * 
     * 
     */
    public List<Sample> getSample() {
        if (sample == null) {
            sample = new ArrayList<Sample>();
        }
        return this.sample;
    }

    /**
     * Gets the value of the name property.
     * 
     * @return
     *     possible object is
     *     {@link String }
     *     
     */
    public String getName() {
        return name;
    }

    /**
     * Sets the value of the name property.
     * 
     * @param value
     *     allowed object is
     *     {@link String }
     *     
     */
    public void setName(String value) {
        this.name = value;
    }

}
