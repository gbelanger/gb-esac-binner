package gb.esac.binner;

import java.util.ArrayList;
import java.util.List;
import java.util.Iterator;
import org.apache.log4j.Logger;


public class DataBinList {

    private static Logger logger  = Logger.getLogger(DataBinList.class);

    private ArrayList<DataBin> dataBinList;

    public DataBinList() {
	dataBinList = new ArrayList<DataBin>();
    }

    public Iterator<DataBin> getIterator() {
	return dataBinList.iterator();
    }


}