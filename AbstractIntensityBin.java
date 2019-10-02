package gb.esac.binner;

import org.apache.log4j.Logger;

/**

The class <code>AbstractIntensityBin</code> extends <code>AbstractBin</code> but implements
 both the IBin and IIntensity interfaces in order to combine the information held by 
the notion of a bin on the horizontal axis and that of an intensity on the vertical 
axis. An <code>AbstractIntensityBin</code> object therefore carries all of the information 
needed to represent a data bin which includes the bin info and the intensity info.

Note that because we can implement several interfaces, and with this inherit all their 
methods, but can only extend one class, if some or all of the interface  methods are 
implemented in abstract classes, only those from the abstrac class that is extended 
will be inherent as implemented there. The methods implemented in any other abstract 
class that implements any other interface will have to be implemented here again. 

This is why we have to implement the methods already implemented in 
<code>AbstractIntensity</code> that relate to setting and getting the intensity and 
the error on it if there one. It's not very elegant, but it's the only way to do it.

 @author <a href="mailto: guilaume.belanger@esa.int">Guillaume Belanger</a>, ESA/ESAC.
 @created March 2013
 @version August 2018

 **/


public abstract class AbstractIntensityBin extends AbstractBin implements IBin, IIntensity {

    private static Logger logger  = Logger.getLogger(AbstractIntensityBin.class);
    // bin attributes
    private double leftEdge;
    private double rightEdge;
    private double[] edges;
    private double width;
    private double centre;
    // values and errors
    private double value = Double.NaN;
    private double error = Double.NaN;
    private double variance = Double.NaN;
    private boolean errorIsSet = false;
    // units and description
    private static String units = null;
    private static String description = null;
    private static boolean unitsAreSet = false;
    private static boolean descriptionIsSet = false;

    //  Methods from IBin are the following:
    //
    // public double[] getEdges();
    // public double getLeftEdge();
    // public double getRightEdge();
    // public double getWidth();
    // public double getCentre();
    // public boolean contains(double value);
    // public boolean contains(IBin bin);
    // public boolean overlaps(IBin bin);
    //
    // We get all of them from extending AbstractBin


    // Constructors are package-private
    protected AbstractIntensityBin() {}

    protected AbstractIntensityBin(AbstractIntensityBin intensityBin) throws BinningException {
	setEdges(intensityBin.getEdges());
	setValue(intensityBin.getValue());
	if (intensityBin.errorIsSet()) {
	    setError(intensityBin.getError());
	}
	printInfo();
    }

    protected AbstractIntensityBin(IIntensity intensity, IBin bin) throws BinningException {
	setEdges(bin.getEdges());
	setValue(intensity.getValue());
	if (intensity.errorIsSet()) {
	    setError(intensity.getError());
	}
	printInfo();	
    }

    //  With value
    protected AbstractIntensityBin(double leftEdge, double rightEdge, double value) throws BinningException {
	setEdges(leftEdge, rightEdge);
	setValue(value);
	printInfo();
    }
    protected AbstractIntensityBin(double leftEdge, double rightEdge, double value, String units) throws BinningException {
	setEdges(leftEdge, rightEdge);
	setValue(value);
	setUnits(units);
	printInfo();
    }
    protected AbstractIntensityBin(double leftEdge, double rightEdge, double value, String units, String description) throws BinningException {
	setEdges(leftEdge, rightEdge);
	setValue(value);
	setUnits(units);
	setDescription(description);
	printInfo();
    }
    
    //  With value and error
    protected AbstractIntensityBin(double leftEdge, double rightEdge, double value, double error) throws BinningException {
	setEdges(leftEdge, rightEdge);
	setValue(value, error);
	printInfo();
    }
    protected AbstractIntensityBin(double leftEdge, double rightEdge, double value, double error, String units) throws BinningException {
	setEdges(leftEdge, rightEdge);
	setValue(value, error);
	setUnits(units);
	printInfo();
    }
    protected AbstractIntensityBin(double leftEdge, double rightEdge, double value, double error, String units, String description) throws BinningException {
	setEdges(leftEdge, rightEdge);
	setValue(value, error);
	setUnits(units);
	setDescription(description);
	printInfo();
    }

    //  Print info
    private void printInfo() {
	logger.info("New IntensityBin is ready: ["+this.getLeftEdge()+", "+this.getRightEdge()+"]");
	logger.info("  Centre = "+this.getCentre());
	logger.info("  Width = "+this.getWidth());
	logger.info("  Intensity = "+this.getValue());
	logger.info("  Error = "+this.getError());
    }
    

    //  Methods duplicated from AbstractIntensity
    
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
    ////  Value and error
    // public IIntensity getIntensity() {

    // }

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

    public boolean errorIsSet() {
	return this.errorIsSet;
    }
    ////  Units and description
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
    
    public boolean unitsAreSet() {
	return this.unitsAreSet;
    }

    public boolean descriptionIsSet() {
	return this.descriptionIsSet;
    }

    
    // Abstract methods
    //public abstract AbstractIntensityBin[] split(double whereToSplit, AbstractIntensityBin previousBin, AbstractIntensityBin nextBin) throws BinningException;
    //public abstract AbstractIntensityBin joinWith(AbstractIntensityBin bin) throws BinningException;

}
