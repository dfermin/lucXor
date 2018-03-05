package lucxor;

/**
 * @author Dmitry Avtonomov
 */
public class Mod {

  public final Site site;
  public final double massDiff;
  public final double mass;
  public final boolean isVariable;

  final boolean isTerminal;
  final boolean isProteinTerm;
  final boolean isPepTerm;

  public Mod(Site site, double massDiff, double mass, boolean isVariable) {
    this.site = site;
    this.massDiff = massDiff;
    this.mass = mass;
    this.isVariable = isVariable;

    isProteinTerm = site == Site.PROT_C_TERM || site == Site.PROT_N_TERM;
    isPepTerm = site == Site.PEP_C_TERM || site == Site.PEP_N_TERM;
    isTerminal = isPepTerm || isProteinTerm;
  }

  public Site getSite() {
    return site;
  }

  public double getMassDiff() {
    return massDiff;
  }

  public double getMass() {
    return mass;
  }

  public boolean isVariable() {
    return isVariable;
  }

  public boolean isTerminal() {
    return isTerminal;
  }

  public boolean isProteinTerm() {
    return isProteinTerm;
  }

  public boolean isPepTerm() {
    return isPepTerm;
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

    if (Double.compare(mod.massDiff, massDiff) != 0) {
      return false;
    }
    if (Double.compare(mod.mass, mass) != 0) {
      return false;
    }
    if (isVariable != mod.isVariable) {
      return false;
    }
    return site == mod.site;
  }

  @Override
  public int hashCode() {
    int result;
    long temp;
    result = site.hashCode();
    temp = Double.doubleToLongBits(massDiff);
    result = 31 * result + (int) (temp ^ (temp >>> 32));
    temp = Double.doubleToLongBits(mass);
    result = 31 * result + (int) (temp ^ (temp >>> 32));
    result = 31 * result + (isVariable ? 1 : 0);
    return result;
  }
}
