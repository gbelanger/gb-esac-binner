package gb.esac.binner;

import org.apache.log4j.Logger;

/**

The class <code>BinTester</code> is written to test <code>Bin</code>.

 @author <a href="mailto: guilaume.belanger@esa.int">Guillaume Belanger</a>, ESA/ESAC.
 @created August 2018
 @version August 2018

 **/

public class BinTester {

    private static Logger logger  = Logger.getLogger(BinTester.class);

    public static void main(String[] args) throws BinningException {

	double number = Double.NaN + 2;
	System.out.println(number);
	System.exit(0);
	

	logger.info("Bin bin1 = new Bin(0,1);");
	Bin bin1 = new Bin(0,1);
	printProperties(bin1);	

	logger.info("bin = new Bin(Double.NaN,Double.NaN);");	
	Bin bin2 = new Bin(Double.NaN,Double.NaN);
	printProperties(bin2);

	logger.info("Bin bin3 = new Bin(bin1);");		
	Bin bin3 = new Bin(bin1);	
	printProperties(bin3);

	logger.info("Bin bin4 = new Bin(2,3);");			
	Bin bin4 = new Bin(2,3);
	printProperties(bin4);

	logger.info("bin4.contains(bin3) = "+bin4.contains(bin3));
	logger.info("bin4.contains(bin1) = "+bin4.contains(bin1));
	logger.info("");
	
	logger.info("Bin bin5 = new Bin(0,2);");
	Bin bin5 = new Bin(0,2);
	printProperties(bin5);
	
	logger.info("bin5.contains(bin1) = "+bin5.contains(bin1));
	logger.info("bin5.contains(bin4) = "+bin5.contains(bin4));
	logger.info("");

	logger.info("Bin bin6 = bin1.joinWith(bin4);");
	//Bin bin6 = bin1.joinWith(bin4);
	Bin bin6 = bin1.joinWith(bin5);
	printProperties(bin6);
	
	logger.info("bin1.overlaps(bin5) = "+bin1.overlaps(bin5));
	logger.info("");

	logger.info("Bin bin7 = new Bin(1,3);");
	Bin bin7 = new Bin(1,3);
	printProperties(bin7);
	logger.info("Bin bin8 = bin1.joinWith(bin7);");	
	Bin bin8 = bin1.joinWith(bin7);
	printProperties(bin8);

	logger.info("Bin[] splitBins = bin1.split(0.5);");
	Bin[] splitBins = bin1.split(0.5);
	printProperties(splitBins[0]);
	printProperties(splitBins[1]);
	
	logger.info("Bin bin9 = splitBins[0].joinWith(splitBins[1]);");
	Bin bin9 = splitBins[0].joinWith(splitBins[1]);	
	printProperties(bin9);
	
    }

    private static void printProperties(Bin bin) {
	logger.info("printProperties(bin):");
	logger.info("  bin.getLeftEdge() = "+bin.getLeftEdge());
	logger.info("  bin.getRightEdge() = "+bin.getRightEdge());
	logger.info("  bin.getCentre() = "+bin.getCentre());
	logger.info("  bin.getWidth() = "+bin.getWidth());
	logger.info("");
    }
}
