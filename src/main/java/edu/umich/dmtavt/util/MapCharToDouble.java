package edu.umich.dmtavt.util;

/**
 * @author Dmitry Avtonomov
 */
public class MapCharToDouble {
  protected double[] map;

  public double get(char ch) {
    return map[(int)ch];
  }

  public double get(int i) {
    return map[i];
  }

  public int size() {
    return map.length;
  }
}
