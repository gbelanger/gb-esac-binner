package gb.esac.binner;


public class BinningUtilsTest {



    public static void main(String[] args) throws BinningException {
	

	double xmin = 0;
	double xmax = 100;
	double binWidth = 10;

	System.out.println("Test 1");
	double[] binEdges = BinningUtils.getBinEdges(xmin, xmax, binWidth);
	for ( int i=0; i < binEdges.length/2; i++ ) {
	    System.out.println(binEdges[2*i] +"\t"+ binEdges[2*i+1]+"\t binWidth = "+(binEdges[2*i+1]-binEdges[2*i]) );
	}


	System.out.println("\nTest 2");	
 	xmin = 100;
 	xmax = 110;
 	binEdges = BinningUtils.getBinEdges(xmin, xmax, binWidth);
	for ( int i=0; i < binEdges.length/2; i++ ) {
	    System.out.println(binEdges[2*i] +"\t"+ binEdges[2*i+1]+"\t binWidth = "+(binEdges[2*i+1]-binEdges[2*i]) );
	}


	System.out.println("\nTest 3");	
	xmin = 0;
	xmax = 1;
	binWidth = 0.1;
	binEdges = BinningUtils.getBinEdges(xmin, xmax, binWidth);
	for ( int i=0; i < binEdges.length/2; i++ ) {
	    System.out.println(binEdges[2*i] +"\t"+ binEdges[2*i+1]+"\t binWidth = "+(binEdges[2*i+1]-binEdges[2*i]) );
	}

	System.out.println("\nTest 4");
	xmin = -Double.MAX_VALUE;
	xmax = Double.MAX_VALUE;
	binWidth = Double.MAX_VALUE/10;
	binEdges = BinningUtils.getBinEdges(xmin, xmax, binWidth);
	for ( int i=0; i < binEdges.length/2; i++ ) {
	    System.out.println(binEdges[2*i] +"\t"+ binEdges[2*i+1]+"\t binWidth = "+(binEdges[2*i+1]-binEdges[2*i]) );
	}


	System.out.println("\nTest 5");
	xmin = 1;
	xmax = 10;
	int nbins = 10;
	binEdges = BinningUtils.getBinEdgesInLogSpace(xmin, xmax, nbins);
	for ( int i=0; i < binEdges.length/2; i++ ) {
	    System.out.println(binEdges[2*i] +"\t"+ binEdges[2*i+1]+"\t binWidth in log space = "+(Math.log(binEdges[2*i+1]) - Math.log(binEdges[2*i])) );
	}


	System.out.println("\nTest 6");
	double[] linearBinCentres = BinningUtils.getBinCentresFromBinEdges(binEdges);
	for ( int i=0; i < binEdges.length/2; i++ ) {
	    System.out.println(binEdges[2*i] +"\t"+ linearBinCentres[i] +"\t"+ binEdges[2*i+1] +"\t binWidth in linear space = "+ (binEdges[2*i+1]-binEdges[2*i])+ "\t distance from edges in linear space = "+(linearBinCentres[i] - binEdges[2*i]) +"\t"+ (binEdges[2*i+1] - linearBinCentres[i]));
	}


    }


}