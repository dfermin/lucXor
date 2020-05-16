package lucxor;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

public class MathFunctionsTest {

  @Test
  public void combinatorial() {

    Assert.assertEquals(MathFunctions.combinatorial(3.345673, 4), 0.25, 0.00001);
  }

  @Test
  public void factorial() {
    Assert.assertEquals(MathFunctions.factorial(4), 24.0, 0.0000000003);
  }

  @Test
  public void roundDouble() {
    Assert.assertEquals(MathFunctions.roundDouble(23.44345, 3), 23.443, 0.000001);
  }

}
