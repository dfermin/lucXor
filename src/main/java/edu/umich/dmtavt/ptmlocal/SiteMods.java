package edu.umich.dmtavt.ptmlocal;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableListMultimap;
import com.google.common.collect.Multimaps;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

/**
 * @author Dmitry Avtonomov
 */
public class SiteMods {
  /** To hold mapping for any of the extended ASCII table symbols. */
  private static final int LEN = Site.MAX_ASCII_INDEX + 1;

  private final Mod[][] map;
  private final boolean isVariable;

  public SiteMods(boolean isVariable) {
    this.map = new Mod[LEN][];
    this.isVariable = isVariable;
  }

  public Mod[] get(char site) {
    return map[(int)site];
  }

  /**
   * @return true if modifications in this mapping are considered fixed.
   */
  public boolean isVariable() {
    return isVariable;
  }

  static class Builder {
    SiteMods sm;
    Set<Mod> mods;

    public Builder(boolean isFixed) {
      sm = new SiteMods(isFixed);
      mods = new HashSet<>();
    }

    public void add(Mod mod) {
      if (mod.isVariable != sm.isVariable)
        throw new IllegalArgumentException("Trying to add a Mod that is variable/fixed to a wrong "
            + "type of mapping SiteMods container builder.");
      mods.add(mod);
    }

    public SiteMods create() {
      Arrays.fill(sm.map, new Mod[0]);
      final ImmutableListMultimap<Site, Mod> index = Multimaps.index(mods, Mod::getSite);
      for (Site site : index.keySet()) {
        final ImmutableList<Mod> list = index.get(site);
        sm.map[site.ch] = list.toArray(new Mod[list.size()]);
      }

      return sm;
    }
  }
}
