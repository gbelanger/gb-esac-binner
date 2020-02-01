package gb.esac.binner;

/** 

The interface <code>IBin</code> defines the methods that any kind of Bin object
must have. A bin is container that holds data, and every bin is defined by its
edges, which define its width and its centre. Therefore, any bin must provide 
methods to access these. A minimal set of methods with explicit names are
defined here. 

These methods are <code>getEdges()</code> that returns both edges,
<code>getLeftEdge()</code> and <code>getRightEdge()</code> that return each 
edge separately, and then <code>getWidth()</code> and <code>getCentre</code>
that are to save the computation of these quantities from the edges.

In addition, to make working with more than one set of bins in combining and 
resampling, for example, a bin object must also provide methods to check if a 
bin contains or overlaps with another bin. These are explicitly named
<code>contains()</code>, two methods that can test for a particular value or for 
an entire bin, and <code>overlaps()</code>.

 @author <a href="mailto: guilaume.belanger@esa.int">Guillaume Belanger</a>, ESA/ESAC.
 @created March 2013
 @version August 2018

**/

public interface IBin {

    public double[] getEdges();
    public double getLeftEdge();
    public double getRightEdge();
    public double getWidth();
    public double getCentre();
    public boolean contains(double value);
    public boolean contains(IBin bin);
    public boolean overlaps(IBin bin);
    public boolean overlapsExactly(IBin bin);
    public boolean isWiderThan(IBin bin);
    public boolean isOfEqualWidth(IBin bin);

    //public IBin[] split(double whereToSplit) throws BinningException;
    //public IBin joinWith(IBin bin) throws BinningException;

}
