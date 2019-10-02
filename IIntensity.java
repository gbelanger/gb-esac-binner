package gb.esac.binner;

/**

The interface <code>IItensity</code> defines the methods that an <code>Intensity</code>
object must provide. These are implemented in <code>AbstractIntensity</code> because
all intensity objects, no matter what kind of intensity they represent, must allow
access to the intensity (<code>getValue()</code>) and to the error and variance
(<code>getError()</code> and <code>getVariance()</code>) on this intensity if it is 
defined (<code>errorIsSet()</code>). There are also two methods to access the units
(<code>getUnits()</code> and a description of the quantity 
(<code>getDescription()</code>).

 @author <a href="mailto: guilaume.belanger@esa.int">Guillaume Belanger</a>, ESA/ESAC.
 @created March 2013
 @version August 2018

 **/

public interface IIntensity {

    double getValue();    
    double getError();
    double getVariance();
    boolean errorIsSet();
    String getUnits();
    String getDescription();    

}
