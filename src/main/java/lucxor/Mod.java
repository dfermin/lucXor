package lucxor;

/**
 * @author Dmitry Avtonomov
 */
public class Mod {

  final Terminus terminus;
  final boolean isProteinEnd;
  final boolean isVariable;
  final double massDiff;

  public Mod(Terminus terminus, boolean isProteinEnd, boolean isVariable, double massDiff) {
    this.terminus = terminus;
    this.isProteinEnd = isProteinEnd;
    this.isVariable = isVariable;
    this.massDiff = massDiff;
  }

  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (o == null || getClass() != o.getClass()) {
      return false;
    }

    Mod mod = (Mod) o;

    if (isProteinEnd != mod.isProteinEnd) {
      return false;
    }
    if (isVariable != mod.isVariable) {
      return false;
    }
    if (Double.compare(mod.massDiff, massDiff) != 0) {
      return false;
    }
    return terminus == mod.terminus;
  }

  @Override
  public int hashCode() {
    int result;
    long temp;
    result = terminus != null ? terminus.hashCode() : 0;
    result = 31 * result + (isProteinEnd ? 1 : 0);
    result = 31 * result + (isVariable ? 1 : 0);
    temp = Double.doubleToLongBits(massDiff);
    result = 31 * result + (int) (temp ^ (temp >>> 32));
    return result;
  }
}
