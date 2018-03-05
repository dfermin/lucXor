package lucxor;

/**
 * @author Dmitry Avtonomov
 */
public class MathHelper {

  private MathHelper() {}

  public static double roundDouble(double value, int numPlaces) {
    double ret;
    double N = Math.pow(10, numPlaces);

    ret = (double) Math.round((value * N)) / N;

    return ret;
  }
}
