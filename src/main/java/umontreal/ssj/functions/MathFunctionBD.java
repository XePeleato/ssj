package umontreal.ssj.functions;

import java.math.BigDecimal;

/**
 * This interface should be implemented by classes which represent univariate
 * mathematical functions. It is used to pass an arbitrary function of one
 * variable as argument to another function. For example, it is used in
 * @ref umontreal.ssj.util.RootFinder to find the zeros of a function.
 *
 * <div class="SSJ-bigskip"></div><div class="SSJ-bigskip"></div>
 */
public interface MathFunctionBD {

    /**
     * Returns the value of the function evaluated at @f$x@f$.
     *  @param x            value at which the function is evaluated
     *  @return function evaluated at `x`
     */
    BigDecimal evaluate (BigDecimal x);

}