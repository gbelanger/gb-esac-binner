package gb.esac.binner;

import org.apache.log4j.Logger;

/**

The abstract class <code>AbstractIntensity</code> implements the methods
defined in the interface <code>IIntensity</code> which allow setting the 
intensity (<code>setIntensity(double intensity)</code>) and the error on this
value (<code>setError(double error)</code>) or both at the same time 
(<code>setIntensity(double intensity, double error)</code>), and then to retrieve 
these (<code>getIntensity()</code>, <code>getError()</code>, and 
<code>getVariance()</code>), but also to check if an error on the intensity 
is set (<code>isErrorSet()</code>). In addition there are also those methods
relating to units (<code>setUnits()</code> and <code>getUnits()</code>) and 
description (<code>setDescription()</code> and <code>getDescription</code>).

 @author <a href="mailto: guilaume.belanger@esa.int">Guillaume Belanger</a>, ESA/ESAC.
 @created March 2013
 @version September 2018

 **/

public abstract class AbstractIntensity implements IIntensity {

    private static Logger logger  = Logger.getLogger(AbstractIntensity.class);
    // values and errors
    private double value = Double.NaN;
    private double error = Double.NaN;
    private double variance = Double.NaN;
    private boolean errorIsSet = false;
    // units and description 
    private String units = null;
    private String description = null;
    private boolean unitsAreSet = false;
    private boolean descriptionIsSet = false;

    
    //  Constructors
    protected AbstractIntensity() {}
    
    protected AbstractIntensity(AbstractIntensity intensity) {
	if (intensity.errorIsSet()) {
	    setValue(intensity.getValue(), intensity.getError());
	}
	else {
	    setValue(intensity.getValue());
	}
	setUnits(intensity.getUnits());
	setDescription(intensity.getDescription());
    }

    protected AbstractIntensity(double value) {
	setValue(value);
    }
    protected AbstractIntensity(double value, String units) {
	setValue(value);
	setUnits(units);
    }
    protected AbstractIntensity(double value, String units, String description) {
	setValue(value);
	setUnits(units);
	setDescription(description);
    }

    protected AbstractIntensity(double value, double error) {
	setValue(value, error);
    }
    protected AbstractIntensity(double value, double error, String units) {
	setValue(value, error);
	setUnits(units);
    }
    protected AbstractIntensity(double value, double error, String units, String description) {
	setValue(value, error);
	setUnits(units);
	setDescription(description);
    }


    //  Protected setters
    private void setValue(double value) {
	this.value = value;
    }

    private void setError(double error) {
	this.error = error;
	if ( !Double.isNaN(error) ) {
	    this.errorIsSet = true;
	    this.variance = Math.pow(error,2);
	}
    }

    private void setValue(double value, double error) {
	setValue(value);
	setError(error);
    }

    
    //  Public setters
    public void setUnits(String units) {
	if (units != null) {
	    this.units = new String(units);
	    this.unitsAreSet = true;
	}
    }

    public void setDescription(String description) {
	if (description != null) {
	    this.description = new String(description);
	    this.descriptionIsSet = true;
	}
    }
    
    //  Public getters    
    public double getValue() {
	return this.value;
    }

    public double getError() {
	if ( !this.errorIsSet() ) {
	    logger.warn("Error is undefined: Returning Double.NaN");
	}
	return this.error;
    }

    public double getVariance() {
	if ( !this.errorIsSet() ) {
	    logger.warn("Variance is undefined: Returning Double.NaN");
	}
	return this.variance;
    }

    public String getUnits() {
	String units = null;
	if (this.unitsAreSet) {
	    units = new String(this.units);
	}
	else {
	    logger.warn("Units are undefined: Returning null String");
	}
	return units;
    }

    public String getDescription() {
	String desc = null;
	if (this.descriptionIsSet) {
	    desc = new String(this.description);
	}
	else {
	    logger.warn("Description is empty: Returning null String");
	}
	return desc;
    }
    
    public boolean errorIsSet() {
	return this.errorIsSet;
    }

    public boolean unitsAreSet() {
	return this.unitsAreSet;
    }

    public boolean descriptionIsSet() {
	return this.descriptionIsSet;
    }

    public void printProperties() {
	logger.info("Properties of Intensity object:");
	logger.info("  Class = "+this.getClass().getSimpleName());
	logger.info("  Value = "+this.getValue());
	if (this.errorIsSet) {
	    logger.info("  Error = "+this.getError());
	    logger.info("  Variance = "+this.getVariance());
	}
	if (this.unitsAreSet) {
	    logger.info("  Units = "+this.getUnits());
	}
	else {
	    logger.info("  Units are undefined");
	}
	if (this.descriptionIsSet) {
	    logger.info("  Description = "+this.getDescription());	
	}
	else {
	    logger.info("  Description is empty");
	}
	logger.info("");
    }


    
}
