package umontreal.ssj.functions;

import java.math.BigDecimal;

public interface MathFunctionBigDec extends MathFunction {
    /**
     * Returns the value of the function evaluated at @f$x@f$.
     *  @param x            value at which the function is evaluated
     *  @return function evaluated at `x`
     */

    BigDecimal evaluate(BigDecimal x);


}
