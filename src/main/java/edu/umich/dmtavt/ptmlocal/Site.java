package edu.umich.dmtavt.ptmlocal;

/**
 * The site in a protein or peptide sequence.
 * @author Dmitry Avtonomov
 */
public enum Site {
  A ('A'),
  B ('B', true),
  C ('C'),
  D ('D'),
  E ('E'),
  F ('F'),
  G ('G'),
  H ('H'),
  I ('I'),
  J ('J', true),
  K ('K'),
  L ('L'),
  M ('M'),
  N ('N'),
  O ('O', true),
  P ('P'),
  Q ('Q'),
  R ('R'),
  S ('S'),
  T ('T'),
  U ('U', true),
  V ('V'),
  W ('W'),
  X ('X', true),
  Y ('Y'),
  Z ('Z', true),
  PEP_N_TERM('<'),
  PEP_C_TERM('>'),
  PROT_N_TERM('['),
  PROT_C_TERM(']');

  char ch;
  boolean isAmbiguous;
  public static final int MAX_ASCII_INDEX = 255;

  Site(char ch, boolean isAmbiguous) {
    this.ch = ch;
    this.isAmbiguous = isAmbiguous;
  }

  Site(char ch) {
    this(ch, false);
  }
}
