package edu.umich.dmtavt.util;

/**
 * @author Dmitry Avtonomov
 */
public class MapCharToChar {
  protected final char[] map;

  public MapCharToChar(char[] map) {
    this.map = map;
  }

  public char get(char ch) {
    return map[(int)ch];
  }

  public char get(int i) {
    return map[i];
  }

  public int size() {
    return map.length;
  }
}
