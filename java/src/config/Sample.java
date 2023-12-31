//
// This file was generated by the JavaTM Architecture for XML Binding(JAXB) Reference Implementation, v2.2.8-b130911.1802 
// See <a href="http://java.sun.com/xml/jaxb">http://java.sun.com/xml/jaxb</a> 
// Any modifications to this file will be lost upon recompilation of the source schema. 
// Generated on: 2020.04.02 at 08:24:17 PM PDT 
//


package config;

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
 *         &lt;element ref="{}Protein"/>
 *         &lt;element ref="{}Concentration" minOccurs="0"/>
 *         &lt;element ref="{}VariableRegionLength"/>
 *         &lt;element ref="{}LeftFlank" minOccurs="0"/>
 *         &lt;element ref="{}RightFlank" minOccurs="0"/>
 *         &lt;element ref="{}LeftBarcode"/>
 *         &lt;element ref="{}RightBarcode"/>
 *         &lt;element ref="{}Round0" minOccurs="0"/>
 *         &lt;element ref="{}Notes"/>
 *       &lt;/sequence>
 *       &lt;attribute name="name" use="required" type="{http://www.w3.org/2001/XMLSchema}string" />
 *       &lt;attribute name="round" use="required" type="{http://www.w3.org/2001/XMLSchema}int" />
 *     &lt;/restriction>
 *   &lt;/complexContent>
 * &lt;/complexType>
 * </pre>
 * 
 * 
 */
@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "", propOrder = {
    "protein",
    "concentration",
    "variableRegionLength",
    "leftFlank",
    "rightFlank",
    "leftBarcode",
    "rightBarcode",
    "round0",
    "notes"
})
@XmlRootElement(name = "Sample")
public class Sample {

    @XmlElement(name = "Protein", required = true, defaultValue = "")
    protected String protein;
    @XmlElement(name = "Concentration", defaultValue = "")
    protected String concentration;
    @XmlElement(name = "VariableRegionLength")
    protected int variableRegionLength;
    @XmlElement(name = "LeftFlank")
    protected String leftFlank;
    @XmlElement(name = "RightFlank")
    protected String rightFlank;
    @XmlElement(name = "LeftBarcode", required = true)
    protected String leftBarcode;
    @XmlElement(name = "RightBarcode", required = true)
    protected String rightBarcode;
    @XmlElement(name = "Round0")
    protected Round0 round0;
    @XmlElement(name = "Notes", required = true, defaultValue = "")
    protected String notes;
    @XmlAttribute(name = "name", required = true)
    protected String name;
    @XmlAttribute(name = "round", required = true)
    protected int round;

    /**
     * Gets the value of the protein property.
     * 
     * @return
     *     possible object is
     *     {@link String }
     *     
     */
    public String getProtein() {
        return protein;
    }

    /**
     * Sets the value of the protein property.
     * 
     * @param value
     *     allowed object is
     *     {@link String }
     *     
     */
    public void setProtein(String value) {
        this.protein = value;
    }

    /**
     * Gets the value of the concentration property.
     * 
     * @return
     *     possible object is
     *     {@link String }
     *     
     */
    public String getConcentration() {
        return concentration;
    }

    /**
     * Sets the value of the concentration property.
     * 
     * @param value
     *     allowed object is
     *     {@link String }
     *     
     */
    public void setConcentration(String value) {
        this.concentration = value;
    }

    /**
     * Gets the value of the variableRegionLength property.
     * 
     */
    public int getVariableRegionLength() {
        return variableRegionLength;
    }

    /**
     * Sets the value of the variableRegionLength property.
     * 
     */
    public void setVariableRegionLength(int value) {
        this.variableRegionLength = value;
    }

    /**
     * Gets the value of the leftFlank property.
     * 
     * @return
     *     possible object is
     *     {@link String }
     *     
     */
    public String getLeftFlank() {
        return leftFlank;
    }

    /**
     * Sets the value of the leftFlank property.
     * 
     * @param value
     *     allowed object is
     *     {@link String }
     *     
     */
    public void setLeftFlank(String value) {
        this.leftFlank = value;
    }

    /**
     * Gets the value of the rightFlank property.
     * 
     * @return
     *     possible object is
     *     {@link String }
     *     
     */
    public String getRightFlank() {
        return rightFlank;
    }

    /**
     * Sets the value of the rightFlank property.
     * 
     * @param value
     *     allowed object is
     *     {@link String }
     *     
     */
    public void setRightFlank(String value) {
        this.rightFlank = value;
    }

    /**
     * Gets the value of the leftBarcode property.
     * 
     * @return
     *     possible object is
     *     {@link String }
     *     
     */
    public String getLeftBarcode() {
        return leftBarcode;
    }

    /**
     * Sets the value of the leftBarcode property.
     * 
     * @param value
     *     allowed object is
     *     {@link String }
     *     
     */
    public void setLeftBarcode(String value) {
        this.leftBarcode = value;
    }

    /**
     * Gets the value of the rightBarcode property.
     * 
     * @return
     *     possible object is
     *     {@link String }
     *     
     */
    public String getRightBarcode() {
        return rightBarcode;
    }

    /**
     * Sets the value of the rightBarcode property.
     * 
     * @param value
     *     allowed object is
     *     {@link String }
     *     
     */
    public void setRightBarcode(String value) {
        this.rightBarcode = value;
    }

    /**
     * Gets the value of the round0 property.
     * 
     * @return
     *     possible object is
     *     {@link Round0 }
     *     
     */
    public Round0 getRound0() {
        return round0;
    }

    /**
     * Sets the value of the round0 property.
     * 
     * @param value
     *     allowed object is
     *     {@link Round0 }
     *     
     */
    public void setRound0(Round0 value) {
        this.round0 = value;
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

    /**
     * Gets the value of the round property.
     * 
     */
    public int getRound() {
        return round;
    }

    /**
     * Sets the value of the round property.
     * 
     */
    public void setRound(int value) {
        this.round = value;
    }

}
