package umontreal.ssj.util;

import java.math.BigDecimal;
import java.math.MathContext;

import static java.math.BigDecimal.ROUND_HALF_UP;
import static java.math.BigDecimal.TEN;

public class BigDecimalUtils {

    private static final int MACH_SCALE = Math.abs(getExponent(BigDecimal.valueOf(getMachineEpsilon()))) + 1;
    public static final BigDecimal TWO = new BigDecimal("2.0");
    public static final BigDecimal THREE = new BigDecimal("3.0");

    /**
     * Compare two {@code BigDecimal}s up to a precision.
     * In other words, if the absolute difference between the two numbers falls below a threshold, they are considered equal.
     *
     * @param n1 a {@code BigDecimal}
     * @param n2 a {@code BigDecimal}
     * @param p  the threshold is <i>1e-p</i>
     * @return -1, 0, or 1 when {@code n1} is numerically less than, equal to, or greater than {@code n2}, respectively
     */
    public static int compare(BigDecimal n1, BigDecimal n2, int p) {
        BigDecimal epsilon = BigDecimal.ONE.movePointLeft(p);
        BigDecimal absDelta = n1.subtract(n2).abs();
        int result = absDelta.compareTo(epsilon);
        return result <= 0 ? 0 : n1.compareTo(n2);
    }

    /**
     * Check if two {@code BigDecimal}s are equal up to a precision.
     *
     * @param n1        a {@code BigDecimal}
     * @param n2        a {@code BigDecimal}
     * @param precision the threshold is <i>1e-p</i>
     * @return {@code true} if the numbers are equal up to a precision
     */
    public static boolean equals(BigDecimal n1, BigDecimal n2, int precision) {
        return compare(n1, n2, precision) == 0;
    }

    /**
     * Compute <i>a</i> to the power of <i>b</i>.
     *
     * @param a a base
     * @param b an exponent
     * @return <i>a<sup>b</sup></i>
     */
    public static BigDecimal pow(BigDecimal a, BigDecimal b) {
        return pow(a, b, MACH_SCALE);
    }

    /**
     * Compute <i>a</i> to the power of <i>b</i>.
     *
     * @param a     a base
     * @param b     an exponent
     * @param scale a precision parameter as in {@link BigDecimal}
     * @return <i>a<sup>b</sup></i>
     */
    public static BigDecimal pow(BigDecimal a, BigDecimal b, int scale) {
        BigDecimal result = b;
        result = result.multiply(log(a, scale));
        result = exp(result, scale);

        return result;
    }

    /**
     * Compute <i>a</i> to the power of <i>n</i>, where <i>n</i> is an integer.
     *
     * @param a a base
     * @param n an integer exponent
     * @return <i>a<sup>n</sup></i>
     */
    public static BigDecimal pow(BigDecimal a, int n) {
        return pow(a, n, MACH_SCALE);
    }

    /**
     * Compute <i>a</i> to the power of <i>n</i>, where <i>n</i> is an integer.
     * This is simply a wrapper around {@link BigDecimal#pow(int)} but handles also negative exponents.
     * Use {@link BigDecimal#pow(int)} for arbitrary precision if the exponent is positive.
     *
     * @param a     a base
     * @param n     an exponent
     * @param scale a precision parameter as in {@link BigDecimal}
     * @return <i>a<sup>n</sup></i>
     */
    public static BigDecimal pow(BigDecimal a, int n, int scale) {
        if (n < 0) {
            return BigDecimal.ONE.divide(pow(a, -n, scale), scale, BigDecimal.ROUND_HALF_EVEN);
        }

        return a.pow(n).setScale(scale, BigDecimal.ROUND_HALF_EVEN);
    }
    //</editor-fold>

    //<editor-fold defaultstate="collapsed" desc="log">
    /**
     * Compute <i>log(x)</i>.
     * The base is <i>e</i>, hence the natural log.
     *
     * @param x a number
     * @return <i>log(x)</i>
     * @see <a href="http://en.wikipedia.org/wiki/Natural_logarithm">Wikipedia: Natural logarithm</a>
     */
    public static BigDecimal log(BigDecimal x) {
        return log(x, MACH_SCALE);
    }

    /**
     * Compute <i>log(x)</i> up to a scale.
     * The base is <i>e</i>, hence the natural log.
     *
     * @param x     a number
     * @param scale a precision parameter as in {@link BigDecimal}
     * @return <i>log(x)</i>
     * @see <a href="http://en.wikipedia.org/wiki/Natural_logarithm">Wikipedia: Natural logarithm</a>
     */
    public static BigDecimal log(BigDecimal x, int scale) {
        BigDecimal logx;

        if (x.compareTo(TEN) < 0) {
            logx = logByNewton(x, scale);
        } else {
            int magnitude = x.precision() - x.scale() - 1;
            BigDecimal m = BigDecimal.valueOf(10).pow(magnitude);//x = 10^m * a
            BigDecimal a = x.divide(m);

            BigDecimal log10 = logByNewton(BigDecimal.valueOf(10), scale);
            BigDecimal loga = logByNewton(a, scale);
            logx = BigDecimal.valueOf(magnitude).multiply(log10).add(loga);
        }

        return logx.setScale(scale, BigDecimal.ROUND_HALF_EVEN);
    }

    /**
     * Compute <i>log(x)</i> up to a scale.
     * The base is <i>e</i>, hence the natural log.
     * <p/>
     * The algorithm solves for {@code y}, by Newton's method,
     * <i>e<sup>y</sup> - x = 0</i>.
     * It works best for "small" <i>x</i>.
     *
     * @param x     <i>x</i>
     * @param scale the scale for the BigDecimal result; a precision parameter
     * @return <i>log(x)</i>
     * @see
     * <ul>
     * <li><a href="http://en.wikipedia.org/wiki/Natural_logarithm">Wikipedia: Natural logarithm</a>
     * <li><a href="http://en.wikipedia.org/wiki/Newton%27s_method">Wikipedia: Newton's method</a>
     * </ul>
     */
    private static BigDecimal logByNewton(BigDecimal x, int scale) {
        if (x.signum() <= 0) {
            throw new IllegalArgumentException("x must be > 0");
        }

        int precision = 2 * scale;//TODO: what is a correct scale for log?
        BigDecimal y = x;
        //a list of initial guesses for known values to speed up the convergence
        if (y.compareTo(BigDecimal.valueOf(10)) == 0) {
            y = BigDecimal.valueOf(2.302585092994);
        } else if (y.compareTo(BigDecimal.valueOf(1)) == 0) {
            y = BigDecimal.ZERO;
        }

        for (BigDecimal delta = BigDecimal.ONE; !equals(delta, BigDecimal.ZERO, scale);) {
            BigDecimal ey = exp(y, precision);
            delta = BigDecimal.ONE.subtract(x.divide(ey, precision, BigDecimal.ROUND_HALF_EVEN));
            y = y.subtract(delta);//y' = y - (1 - x / e(y))
        }

        return y.setScale(scale, BigDecimal.ROUND_HALF_EVEN);
    }
    //</editor-fold>

    //<editor-fold defaultstate="collapsed" desc="exp">
    /**
     * Compute <i>e<sup>x</sup></i>.
     *
     * @param x the exponent
     * @return <i>e<sup>x</sup></i>
     * @see <a href="http://en.wikipedia.org/wiki/Exponential_function">Wikipedia: Exponential function</a>
     */
    public static BigDecimal exp(double x) {
        return exp(x, MACH_SCALE);
    }

    /**
     * Compute <i>e<sup>x</sup></i>.
     *
     * @param x     the exponent
     * @param scale a precision parameter as in {@link BigDecimal}
     * @return <i>e<sup>x</sup></i>
     * @see <a href="http://en.wikipedia.org/wiki/Exponential_function">Wikipedia: Exponential function</a>
     */
    public static BigDecimal exp(double x, int scale) {
        return exp(BigDecimal.valueOf(x), scale);
    }

    /**
     * Compute <i>e<sup>x</sup></i>.
     *
     * @param x the exponent
     * @return <i>e<sup>x</sup></i>
     * @see <a href="http://en.wikipedia.org/wiki/Exponential_function">Wikipedia: Exponential function</a>
     */
    public static BigDecimal exp(BigDecimal x) {
        return exp(x, MACH_SCALE);
    }

    /**
     * Get the integral part of a number (discarding the fractional part).
     *
     * @param num a {@code BigDecimal}
     * @return the integral part of the number
     */
    public static BigDecimal getWhole(BigDecimal num) {
        BigDecimal result = num.setScale(0, BigDecimal.ROUND_DOWN);
        return result;
    }

    /**
     * Get the fractional part of a number.
     * This is the same as the number subtracting the whole part.
     * For a -ve. number, the fractional part is also -ve.
     * For example, for <i>-3.1415</i>, the whole is <i>-3</i> and the fractional part is <i>-0.1415</i>.
     *
     * @param num a {@code BigDecimal}
     * @return the fractional part of the number
     */
    public static BigDecimal getFractional(BigDecimal num) {
        BigDecimal result = num.subtract(getWhole(num));
        return result;
    }

    /**
     * Compute <i>e<sup>x</sup></i>.
     *
     * @param x     the exponent
     * @param scale a precision parameter as in {@link BigDecimal}
     * @return <i>e<sup>x</sup></i>
     * @see <a href="http://en.wikipedia.org/wiki/Exponential_function">Wikipedia: Exponential function</a>
     */
    public static BigDecimal exp(BigDecimal x, int scale) {
        int whole = getWhole(x).intValue();
        BigDecimal fraction = getFractional(x);//x = whole + fraction

        BigDecimal e = expByTaylor(BigDecimal.ONE, scale);
        BigDecimal exp1 = pow(e, whole, scale);//e ^ whole
        BigDecimal exp2 = expByTaylor(fraction, scale);//e ^ fraction
        BigDecimal result = exp1.multiply(exp2);

        return result.setScale(scale, BigDecimal.ROUND_HALF_EVEN);
    }

    /**
     * Compute <i>e<sup>x</sup></i> using Taylor series expansion.
     * This works best for small <i>x</i> as the convergence is quick.
     *
     * @param x     the exponent
     * @param scale a precision parameter as in {@link BigDecimal}
     * @return <i>e<sup>x</sup></i>
     * @see <a href="http://en.wikipedia.org/wiki/Exponential_function">Wikipedia: Exponential function</a>
     */
    private static BigDecimal expByTaylor(BigDecimal x, int scale) {
        final MathContext mc = new MathContext(scale);
        BigDecimal result = BigDecimal.ZERO, newResult = BigDecimal.ONE;

        BigDecimal term = new BigDecimal(1, mc);
        for (int i = 1; !equals(result, newResult, scale); ++i) {
            result = newResult;

            term = term.multiply(x, mc);//x ^ i
            term = term.divide(BigDecimal.valueOf(i), scale + 1, BigDecimal.ROUND_HALF_EVEN);//(x ^ i) / i!
            newResult = result.add(term);//∑((x ^ i) / i!)
        }

        return result.setScale(scale, BigDecimal.ROUND_HALF_EVEN);
    }

    /**
     * The machine epsilon differs per machine so we need to compute it.
     * <p/>
     * This algorithm does not actually determine the machine epsilon.
     * Rather, it determines a number within a factor of two (one order of magnitude) of the true machine epsilon,
     * using a linear search.
     *
     * @return the machine epsilon
     * @see <a href="http://en.wikipedia.org/wiki/Machine_epsilon#Approximation_using_Java">Wikipedia: Approximation using Java</a>
     */
    private static double getMachineEpsilon() {
        double epsilon = 1;

        do {
            epsilon /= 2;
        } while ((1.0 + (epsilon / 2.0)) != 1.0);

        return epsilon;
    }

    private static int getExponent(BigDecimal x) {
        if (x.equals(BigDecimal.ZERO)) {
            return 0;
        }

        int b = 0;
        BigDecimal absx = x.abs();
        for (; absx.compareTo(TEN) >= 0;) {//|x| >= 10
            ++b;
            absx = absx.divide(TEN);
        }

        for (; absx.compareTo(BigDecimal.ONE) < 0;) {//|x| < 1
            --b;
            absx = absx.multiply(TEN);
        }

        return b;
    }

    public static BigDecimal sqrt(BigDecimal A) {
        BigDecimal x0 = BigDecimal.ZERO;
        BigDecimal x1 = new BigDecimal(Math.sqrt(A.doubleValue()));
        while (!x0.equals(x1)) {
            x0 = x1;
            x1 = A.divide(x0, MACH_SCALE, ROUND_HALF_UP);
            x1 = x1.add(x0);
            x1 = x1.divide(TWO, MACH_SCALE, ROUND_HALF_UP);

        }
        return x1;
    }

    public static BigDecimal nthRoot(final int n, final BigDecimal a) {
        return nthRoot(n, a, BigDecimal.valueOf(.1).movePointLeft(MACH_SCALE));
    }

    // https://stackoverflow.com/questions/22695654/computing-the-nth-root-of-p-using-bigdecimals
    private static BigDecimal nthRoot(final int n, final BigDecimal a, final BigDecimal p) {
        if (a.compareTo(BigDecimal.ZERO) < 0) {
            throw new IllegalArgumentException("nth root can only be calculated for positive numbers");
        }
        if (a.equals(BigDecimal.ZERO)) {
            return BigDecimal.ZERO;
        }
        BigDecimal xPrev = a;
        BigDecimal x = a.divide(new BigDecimal(n), MACH_SCALE, ROUND_HALF_UP);  // starting "guessed" value...
        while (x.subtract(xPrev).abs().compareTo(p) > 0) {
            xPrev = x;
            x = BigDecimal.valueOf(n - 1.0)
                    .multiply(x)
                    .add(a.divide(x.pow(n - 1), MACH_SCALE, ROUND_HALF_UP))
                    .divide(new BigDecimal(n), MACH_SCALE, ROUND_HALF_UP);
        }
        return x;
    }

    public static double[] toDoubleArray(BigDecimal[] bigDecimals) {
        double[] result = new double[bigDecimals.length];
        for (int i = 0; i < bigDecimals.length; i++) {
            result[i] = bigDecimals[i].doubleValue();
        }
        return result;
    }

    public static double[][] toDoubleMatrix(BigDecimal[][] matrix) {
        double[][] ret = new double[matrix.length][];
        for (int i = 0; i < ret.length; i++) {
            ret[i] = new double[matrix[i].length];
            for (int j = 0; j < matrix[i].length; j++)
                ret[i][j] = matrix[i][j].doubleValue();
        }

        return ret;
    }

    public static BigDecimal[] toDecimalArray(double[] doubles, MathContext mc) {
        BigDecimal[] result = new BigDecimal[doubles.length];
        for (int i = 0; i < doubles.length; i++) {
            result[i] = new BigDecimal(doubles[i], mc);
        }
        return result;
    }

    public static BigDecimal[][] toDecimalMatrix(double[][] matrix) {
        BigDecimal[][] ret = new BigDecimal[matrix.length][];
        for (int i = 0; i < ret.length; i++) {
            ret[i] = new BigDecimal[matrix.length];
            for (int j = 0; j < matrix[i].length; i++)
                ret[i][j] = BigDecimal.valueOf(matrix[i][j]);
        }

        return ret;
    }

}
