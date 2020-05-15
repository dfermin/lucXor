package lucxor;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.*;

public class MathFunctionsTest {

  @Test
  public void combinatorial() {

    Assert.assertEquals(MathFunctions.combinatorial(3.345673, 4), 0.25, 0.00001);
  }
}
