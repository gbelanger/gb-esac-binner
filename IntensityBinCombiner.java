package gb.esac.binner;

import org.apache.log4j.Logger;

public final class IntensityBinCombiner {

    private static Logger logger  = Logger.getLogger(IntensityBinCombiner.class);
    
    // For DensityBin
    
    public static DensityBin join(DensityBin bin1, DensityBin bin2) throws BinningException {
	logger.info("Attempting to join density bins...");
	DensityBin newDensityBin = null;
	double left1 = bin1.getLeftEdge();
	double right1 = bin1.getRightEdge();
	double left2 = bin2.getLeftEdge();
	double right2 = bin2.getRightEdge();
	double leftEdge = Math.min(left1, left2);
	double rightEdge = Math.max(right1, right2);
	// Adjacent bins 
	if (right1 == left2 || left1 == right2) {
	    if (right1 == left2) {
		logger.info("  Bins are adjacent (bin1 is to the left of bin2).");
	    }
	    else {
		logger.info("  Bins are adjacent (bin1 is to the right of bin2).");
	    }
	}
	// Bin1 contains the other
    	else if (bin1.contains(bin2) || bin2.contains(bin1)) {
	    if (bin1.contains(bin2)) {
		logger.info("  Bin1 contains bin2.");
	    }
	    else {
		logger.info("  Bin2 contains bin1.");
	    }
    	}
	// They overlap (partially or completely)
    	else if (bin1.overlaps(bin2)) {
	    logger.info("  Bins overlap.");	
    	}
	// Bins are disjointed
	else {
	    throw new BinningException("Cannot join bins: ["+left1+", "+right1+"] and ["+left2+", "+right2+"] are disjointed.");
	}
	Density density = calculateCombinedIntensity(bin1, bin2);
	newDensityBin = new DensityBin(leftEdge, rightEdge, density.getValue(), density.getError(), density.getUnits(), density.getDescription());
	logger.info("Joining complete");
    	return newDensityBin;
    }

    private static Density calculateCombinedIntensity(DensityBin bin1, DensityBin bin2) throws BinningException {
	if (bin1.getUnits() != bin2.getUnits()) {
	    throw new BinningException("Cannot combine intensities: Units ("+bin1.getUnits()+" and "+bin2.getUnits()+") are not the same.");
	}
	String units = bin1.getUnits();
	String description = "Combined density";
	double newError = Double.NaN;
	double newValue = 0;
	if (bin1.errorIsSet() || bin2.errorIsSet()) {
	    if (bin1.errorIsSet() && bin2.errorIsSet()) {
		newError = getErrorFromWeights(bin1.getVariance(), bin2.getVariance());
		// Use variances as weights
		description = description+" using variances as weights.";
		newValue = combineValues(bin1.getValue(), bin1.getVariance(), bin2.getValue(), bin2.getVariance());
	    }
	    else {
		if (bin1.errorIsSet()) {newError = bin1.getError();}
		else {newError = bin2.getError();}
		// Use bin widths as weights
		description = description+" using bin widths as weights.";
		newValue = combineValues(bin1.getValue(), bin1.getWidth(), bin2.getValue(), bin2.getWidth());
	    }
	}
	else {
	    description = description+" using bin widths as weights.";	    
	    newValue = combineValues(bin1.getValue(), bin1.getWidth(), bin2.getValue(), bin2.getWidth());	    
	}
	//  Add previous descriptions from each bin
	if (bin1.getDescription() != null) {
	    description = description+" "+bin1.getDescription()+" (bin1);";
	    if (bin2.getDescription() != null) {
		description = description+" "+bin2.getDescription()+" (bin2).";
	    }
	}
	return new Density(newValue, newError);
    }


    // For AbsoluteQuantityBin
    
    public static AbsoluteQuantityBin join(AbsoluteQuantityBin bin1, AbsoluteQuantityBin bin2) throws BinningException {
	logger.info("Trying to join absolute quantity bins...");
	AbsoluteQuantityBin newAbsoluteQuantityBin = null;
	double left1 = bin1.getLeftEdge();
	double right1 = bin1.getRightEdge();
	double left2 = bin2.getLeftEdge();
	double right2 = bin2.getRightEdge();
	double leftEdge = Math.min(left1, left2);
	double rightEdge = Math.max(right1, right2);
	// Adjacent bins 
	if (right1 == left2 || left1 == right2) {
	    if (right1 == left2) {
		logger.info("  Bins are adjacent (bin1 is to the left of bin2).");
	    }
	    else {
		logger.info("  Bins are adjacent (bin1 is to the right of bin2).");
	    }
	}
	// Bin1 contains the other
    	else if (bin1.contains(bin2) || bin2.contains(bin1)) {
	    if (bin1.contains(bin2)) {
		logger.info("  Bin1 contains bin2.");
	    }
	    else {
		logger.info("  Bin2 contains bin1.");
	    }
    	}
	// They overlap (partially or completely)
    	else if (bin1.overlaps(bin2)) {
	    logger.info("  Bins overlap.");	
    	}
	// Bins are disjointed
	else {
	    throw new BinningException("Cannot join bins: ["+left1+", "+right1+"] and ["+left2+", "+right2+"] are disjointed.");
	}
	AbsoluteQuantity quantity = calculateCombinedIntensity(bin1, bin2);
	newAbsoluteQuantityBin = new AbsoluteQuantityBin(leftEdge, rightEdge, quantity.getValue());
	logger.info("Joining complete.");
    	return newAbsoluteQuantityBin;
    }

    private static AbsoluteQuantity calculateCombinedIntensity(AbsoluteQuantityBin bin1, AbsoluteQuantityBin bin2) throws BinningException {
	if (bin1.getUnits() != bin2.getUnits()) {
	    throw new BinningException("Cannot combine intensities: Units ("+bin1.getUnits()+" and "+bin2.getUnits()+") are not the same");
	}
	String units = bin1.getUnits();
	String description = "Combined absolute quantity using bin widths as weights.";
	if (bin1.getDescription() != null) {
	    description = description+" "+bin1.getDescription()+" (bin1);";
	    if (bin2.getDescription() != null) {
		description = description+" "+bin2.getDescription()+" (bin2).";
	    }
	}
	double newValue = combineValues(bin1.getValue(), bin1.getWidth(), bin2.getValue(), bin2.getWidth());
	return new AbsoluteQuantity(newValue, units, description);
    }

    
    //  Private

    private static double combineValues(double value1, double weight1, double value2, double weight2) {
	double newDensityValue = Double.NaN;
	if (!Double.isNaN(new Double(value1)) || !Double.isNaN(new Double(value2))) {
	    if (!Double.isNaN(new Double(value1)) && !Double.isNaN(new Double(value2))) {
		double weightedSum = weight1*value1 + weight2*value2;
		double sumOfWeights = weight1 + weight2;
		newDensityValue = weightedSum/sumOfWeights;
	    }
	    else {
		if (Double.isNaN(new Double(value1))) {newDensityValue = value2;}
		else {newDensityValue = value1;}
	    }
	}
	return newDensityValue;
    }

    private static double getErrorFromWeights(double w1, double w2) {
	double newError = 0;
	if (w1 != 0 || w2 != 0) {
	    double sumOfWeights = w1 + w2;
	    newError = 1./Math.sqrt(sumOfWeights);
	}
	return newError;
    }
    

}
