package gb.esac.binner;

import java.util.ArrayList;
import java.util.List;
import java.util.Iterator;
import org.apache.log4j.Logger;


public class IntensityList {

    private static Logger logger  = Logger.getLogger(IntensityList.class);

    private ArrayList<Intensity> intensityList;

    public IntensityList() {
	intensityList = new ArrayList<Intensity>();
    }

    public Iterator<Intensity> getIterator() {
	return intensityList.iterator();
    }

}