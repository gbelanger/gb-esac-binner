package gb.esac.binner;

import org.apache.log4j.Logger;

/**

The <code>Density</code> class represents a measurement of any kind of quantity.
All the methods it must implement are implemented in <code>AbstractIntensity</code>.

 @author <a href="mailto: guilaume.belanger@esa.int">Guillaume Belanger</a>, ESA/ESAC.
 @created March 2013
 @version August 2018

 **/

public class Density extends AbstractIntensity {

    private static Logger logger  = Logger.getLogger(Density.class);
	
    //  Constructors
    private Density() {
	super();
    }

    public Density(Density density) {
	super(density);
    }

    //  with value
    public Density(double value) {
	super(value);
    }
    public Density(double value, String units) {
	super(value, units);
    }
    public Density(double value, String units, String description) {
	super(value, units, description);
    }

    //  with value and error
    public Density(double value, double error) {
	super(value, error);
    }
    public Density(double value, double error, String units) {
	super(value, error, units);
    }
    public Density(double value, double error, String units, String description) {
	super(value, error, units, description);
    }
    

}
