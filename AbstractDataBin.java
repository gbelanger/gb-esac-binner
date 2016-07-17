package gb.esac.binner;

import org.apache.log4j.Logger;


public abstract class AbstractDataBin extends AbstractBin implements IBin, IIntensity  {

    private static Logger logger  = Logger.getLogger(AbstractDataBin.class);

    protected double centre;
    protected double width;
    protected double leftEdge;
    protected double rightEdge;
    protected double[] edges;

    protected double intensity = Double.NaN;
    protected double error = Double.NaN;
    protected double variance = Double.NaN;
    protected boolean errorIsSet = false;


    //  Protected methods for Intensity defined here

    protected void setIntensity(double intensity) {
	this.intensity = intensity;
    }

    protected void setError(double error) {
	this.error = error;
	if ( !Double.isNaN(error) ) {
	    this.errorIsSet = true;
	    this.variance = error*error;
	}
    }


    //  Methods from IIntensity that must be implemented
    
    @Override
    public double getIntensity() {
	return this.intensity;
    }

    @Override
    public double getError() {
	if ( !this.errorIsSet ) logger.warn("Error is not defined: Returning Double.NaN");
	return this.error;
    }

    @Override
    public double getVariance() {
	if ( !this.errorIsSet ) logger.warn("Error is not defined: Returning Double.NaN");
	return this.variance;
    }

    @Override
    public boolean isErrorSet() {
	return this.errorIsSet;
    }


}