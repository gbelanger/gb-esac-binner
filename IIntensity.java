package gb.esac.binner;

public interface IIntensity {

    double getIntensity();
    double getError();
    double getVariance();
    boolean isErrorSet();

}