package gb.esac.binner;

public interface IBin {

    double[] getEdges();
    double getLeftEdge();
    double getRightEdge();
    double getWidth();
    double getCentre();
    boolean covers(double value);
    boolean covers(IBin bin);
    boolean overlaps(IBin bin);

}
