package gb.esac.binner;

import java.util.ArrayList;
import java.util.Iterator;
import org.apache.log4j.Logger;

/**

The class <code>BinList</code> that holds an <code>ArrayList</code> of 
<code>Bin</code> objects allowing using it with any kind of bin.

 @author <a href="mailto: guilaume.belanger@esa.int">Guillaume Belanger</a>, ESA/ESAC.
 @created March 2013
 @version August 2018

 **/

public class BinList {

    private static Logger logger  = Logger.getLogger(BinList.class);

    private ArrayList<IBin> binList;

    public BinList() {
	binList = new ArrayList<IBin>();
    }

    public Iterator<IBin> getIterator() {
	return binList.iterator();
    }

}
