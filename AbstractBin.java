package gb.esac.binner;

import org.apache.log4j.Logger;

/**

The class <code>AbstractBin</code> is the general representation of a bin object.
It implements the IBin interface that defines the minimal set of methods that any 
bin object must provide.

Because all bins must have certain properties related to the definition of their
edges, all the methods defined in the interface <code>IBin</code> for setting and 
getting the bin edges, as well as those to check for coverage and overlap, are 
implemented here.

Note that all is defined and provided here refers only to the horizontal dimension 
of the bin: only to the container, but nothing about the contents, which must be 
defined elsewhere in an object that implements the <code>IIntensity</code> 
interface as does the <code<Intensity</code> object which extends 
<code>AbstractIntensity</code>.

 @author <a href="mailto: guilaume.belanger@esa.int">Guillaume Belanger</a>, ESA/ESAC.
 @created March 2013
 @version August 2018

 **/

public abstract class AbstractBin implements IBin {

    //  Class variable
    private static Logger logger  = Logger.getLogger(AbstractBin.class);

    //  Instance variables
    private double leftEdge;
    private double rightEdge;
    private double[] edges;
    private double width;
    private double centre;

    //  Constructors are Package-private
    AbstractBin() {}

    AbstractBin(IBin bin) throws BinningException {
    	setEdges(bin.getLeftEdge(), bin.getRightEdge());
    	//printInfo();
    }

    AbstractBin(double leftEdge, double rightEdge) throws BinningException {
    	setEdges(leftEdge, rightEdge);
    	//printInfo();
    }

    //  Print info
    private void printInfo() {
    	logger.info("New Bin is ready: ["+this.getLeftEdge()+", "+this.getRightEdge()+"]");
    	logger.info("  Centre = "+this.getCentre());
    	logger.info("  Width = "+this.getWidth());
    }


    //  package-private setters
    void setEdges(double leftEdge, double rightEdge) throws BinningException {
    	this.leftEdge = leftEdge;
    	this.rightEdge = rightEdge;
    	this.edges = new double[] {leftEdge, rightEdge};
    	setCentre((leftEdge+rightEdge)/2d);
    	setWidth(rightEdge - leftEdge);
    }

    void setEdges(double[] edges) throws BinningException {
    	setEdges(edges[0], edges[1]);
    }

    void setCentre(double centre) {
    	this.centre = centre;
    }

    void setWidth(double width) throws BinningException {
    	if (width <= Math.ulp(0.0)) {
    	    throw new BinningException("Cannot construct zero (or negative) size bin");
    	}
    	this.width = width;
    }


    //  Public methods from IBin that must be implemented
    public double getLeftEdge() {
    	return this.leftEdge;
    }

    public double getRightEdge() {
    	return this.rightEdge;
    }

    public double[] getEdges() {
    	return this.edges;
    }

    public double getWidth() {
    	return this.width;
    }

    public double getCentre() {
    	return this.centre;
    }

    public boolean contains(double value) {
    	boolean boundedFromLeft = value > (this.leftEdge - Math.ulp(this.leftEdge));
    	boolean boundedFromRight = value < (this.rightEdge + Math.ulp(this.rightEdge));
    	return (boundedFromLeft && boundedFromRight);
    }

    public boolean contains(IBin bin) {
    	boolean containsLeftEdge = contains(bin.getLeftEdge());
    	boolean containsRightEdge = contains(bin.getRightEdge());
    	return (containsLeftEdge && containsRightEdge);
    }

    public boolean overlaps(IBin bin) {
    	boolean containsLeftEdge = contains(bin.getLeftEdge());
    	boolean containsRightEdge = contains(bin.getRightEdge());
    	return (containsLeftEdge || containsRightEdge);
    }

    public boolean overlapsExactly(IBin bin) {
        return (this.leftEdge == bin.getLeftEdge() && this.rightEdge == bin.getRightEdge());
    }

    public boolean isWiderThan(IBin bin) {
        return (this.width > bin.getWidth());
    }

    public boolean isOfEqualWidth(IBin bin){
        return (this.width == bin.getWidth());
    }


    //public abstract AbstractBin[] split(double whereToSplit) throws BinningException;
    //public abstract AbstractBin joinWith(AbstractBin bin) throws BinningException;
    
}
