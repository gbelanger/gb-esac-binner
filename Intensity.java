package gb.esac.binner;

import org.apache.log4j.Logger;


public class Intensity extends AbstractIntensity {

    private static Logger logger  = Logger.getLogger(Intensity.class);

    //  Constructors

    private Intensity() {

    }

    public Intensity(Intensity intensity) {
	setIntensity(intensity.getIntensity(), intensity.getError());
    }

    public Intensity(double intensity) {
	setIntensity(intensity);
    }

    public Intensity(double intensity, double error) {
	setIntensity(intensity);
	setError(error);
    }

}