package gb.esac.binner;

import java.util.ArrayList;
import java.util.List;
import java.util.Iterator;
import org.apache.log4j.Logger;


public class BinList {

    private static Logger logger  = Logger.getLogger(BinList.class);

    private ArrayList<Bin> binList;

    public BinList() {
	binList = new ArrayList<Bin>();
    }

    public Iterator<Bin> getIterator() {
	return binList.iterator();
    }

    public ArrayList<Bin> getBinListCopy() {
	return new ArrayList<Bin>(binList);
    }

}