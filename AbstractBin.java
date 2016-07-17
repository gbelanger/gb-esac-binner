package gb.esac.binner;

import org.apache.log4j.Logger;


public abstract class AbstractBin implements IBin {

    //  Class variable
    private static Logger logger  = Logger.getLogger(AbstractBin.class);

    //  Instance variables
    double centre;
    double width;
    double leftEdge;
    double rightEdge;
    double[] edges;


    //  Protected set methods

    protected void setEdges(double leftEdge, double rightEdge) {
	this.leftEdge = leftEdge;
	this.rightEdge = rightEdge;
	this.edges = new double[] {leftEdge, rightEdge};
	setCentre((leftEdge+rightEdge)/2d);
	setWidth(rightEdge - leftEdge);
    }

    protected void setEdges(double[] edge) {
	setEdges(edges[0], edges[1]);
    }

    protected void setCentre(double centre) {
	this.centre = centre;
    }

    protected void setWidth(double width) {
	this.width = width;
    }



    //  Public get methods from IBin that must be implemented
    @Override
    public double[] getEdges() {
	return this.edges;
    }

    @Override
    public double getLeftEdge() {
	return this.leftEdge;
    }

    @Override
    public double getRightEdge() {
	return this.rightEdge;
    }

    @Override
    public double getWidth() {
	return this.width;
    }

    @Override
    public double getCentre() {
	return this.centre;
    }


    //  Other public methods
    @Override
    public boolean covers(double value) {

	boolean valueIsGreaterThanLeftEdge = value > (this.leftEdge-Math.ulp(this.leftEdge));
	boolean valueIsLessThanRightEdge = value < (this.rightEdge+Math.ulp(this.rightEdge));
	return (  valueIsGreaterThanLeftEdge && valueIsLessThanRightEdge );
    }

    @Override
    public boolean covers(IBin bin) {
	double left = bin.getLeftEdge();
	double right = bin.getRightEdge();
	boolean containsLeftEdge = this.leftEdge < (left - Math.ulp(left));
	boolean containsRightEdge = this.rightEdge > (right + Math.ulp(right));
	return (containsLeftEdge && containsRightEdge);
    }

    @Override
    public boolean overlaps(IBin bin) {
	double left = bin.getLeftEdge();
	double right = bin.getRightEdge();
	boolean coversLeftEdgeOfOtherBin = this.rightEdge >= (left + Math.ulp(left));
	boolean coversRightEdgeOfOtherBin = this.leftEdge <=  (right - Math.ulp(right));
	return (coversLeftEdgeOfOtherBin || coversRightEdgeOfOtherBin);
    }

}