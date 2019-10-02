package gb.esac.binner;

import org.apache.log4j.Logger;

/**

The <code>AbsoluteQuantity</code> class represents a measurement of any kind of quantity.
All the methods it must implement are implemented in <code>AbstractAbsoluteQuantity</code>.

 @author <a href="mailto: guilaume.belanger@esa.int">Guillaume Belanger</a>, ESA/ESAC.
 @created March 2013
 @version September 2018

 **/

public class AbsoluteQuantity extends AbstractIntensity {

    private static Logger logger  = Logger.getLogger(AbsoluteQuantity.class);
    
    //  Constructors
    private AbsoluteQuantity() {
	super();
    }
	
    public AbsoluteQuantity(AbsoluteQuantity absoluteQuantity) {
	super(absoluteQuantity);
    }

    public AbsoluteQuantity(double value) {
	super(value);
    }
    public AbsoluteQuantity(double value, String units) {
	super(value, units);
    }
    public AbsoluteQuantity(double value, String units, String description) {
	super(value, units, description);
    }

}
