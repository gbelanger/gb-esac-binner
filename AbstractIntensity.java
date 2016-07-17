package gb.esac.binner;

import org.apache.log4j.Logger;


public abstract class AbstractIntensity implements IIntensity {

    private static Logger logger  = Logger.getLogger(AbstractIntensity.class);

    protected double intensity = Double.NaN;
    protected double error = Double.NaN;
    protected double variance = Double.NaN;
    protected boolean errorIsSet = false;


    //  Protected set methods defined here

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


    //  Public methods from IIntensity that must be implemented
    
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
    public boolean isErrorSet() {
	return this.errorIsSet;
    }

}