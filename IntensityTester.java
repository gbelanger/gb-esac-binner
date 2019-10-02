package gb.esac.binner;

import org.apache.log4j.Logger;

/**

The class <code>IntensityTester</code> is written to test <code>Intensity</code>.

 @author <a href="mailto: guilaume.belanger@esa.int">Guillaume Belanger</a>, ESA/ESAC.
 @created August 2018
 @version September 2018

 **/

public class IntensityTester {

    private static Logger logger  = Logger.getLogger(IntensityTester.class);
    
    public static void main(String[] args) {

	// Density
	
	// // with value only
	double value = 10;
	logger.info("Density int1 = new Density(value);");
	Density int1 = new Density(value);
	printProperties(int1);
	int1.printProperties();
	
	value = 10;
	String units = "units";
	logger.info("Density int1 = new Density(value,units);");
	int1 = new Density(value,units);
	printProperties(int1);
	int1.printProperties();	

	value = 10;
	units = "units";
	String desc = "description";
	logger.info("Density int1 = new Density(value,units,description);");
	int1 = new Density(value,units,desc);
	printProperties(int1);
	int1.printProperties();	

	
	// // with value and error
	double error = 2;
	logger.info("Density int2 = new Density(value,error);");	
	Density int2 = new Density(value,error);
	printProperties(int2);
	int2.printProperties();	

	logger.info("Density int3 = new Density(int2);");
	Density int3 = new Density(int2);
	printProperties(int3);
	int3.printProperties();
	
	units = "cps";
	int3 = new Density(value,error,units);
	printProperties(int3);
	int3.printProperties();
	
	units = "cps2";
	desc = "other description";
	int3 = new Density(value,error,units,desc);
	printProperties(int3);
	int3.printProperties();	
	
	// //  Check public setters
	Density int5 = new Density(6, 2);
	printProperties(int5);
	int5.printProperties();
	int5.setUnits("new units");
	int5.setDescription("new description");
	printProperties(int5);
	int5.printProperties();

	//  AbsoluteQuantity
	AbsoluteQuantity int6 = new AbsoluteQuantity(6);
	int6.printProperties();

	int6 = new AbsoluteQuantity(5,"counts");
	int6.printProperties();

	int6 = new AbsoluteQuantity(4,"counts","X-rays 2-10 keV");
	int6.printProperties();

	AbsoluteQuantity int7 = new AbsoluteQuantity(int6);
	int7.printProperties();
	
    }

    private static void printProperties(AbstractIntensity intensity) {
	logger.info("printProperties(intensity ("+intensity.getClass().getSimpleName()+"))");
	logger.info("  intensity.getValue() = "+intensity.getValue());
	logger.info("  intensity.getError() = "+intensity.getError());
	logger.info("  intensity.getVariance() = "+intensity.getVariance());
	logger.info("  intensity.errorIsSet() = "+intensity.errorIsSet());
	logger.info("  intensity.unitsAreSet() = "+intensity.unitsAreSet());
	logger.info("  intensity.getUnits() = "+intensity.getUnits());
	logger.info("  intensity.descriptionIsSet() = "+intensity.descriptionIsSet());
	logger.info("  intensity.getDescription() = "+intensity.getDescription());	
	logger.info("");
    }
    
}
