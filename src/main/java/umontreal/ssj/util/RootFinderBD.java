package umontreal.ssj.util;

import umontreal.ssj.functions.MathFunction;
import umontreal.ssj.functions.MathFunctionBD;

import java.math.BigDecimal;
import java.math.MathContext;

/**
 * This class provides methods to solve non-linear equations.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class RootFinderBD {
    private static final BigDecimal MINVAL = new BigDecimal("5.0e-308");
    /**
     * Difference between 1.0 and the smallest `double` greater than 1.0.
     */
    public static final BigDecimal DBL_EPSILON = new BigDecimal("2.2204460492503131e-16");
    private RootFinderBD() {}

    /**
     * Computes a root @f$x@f$ of the function in `f` using the
     * Brent-Dekker method. The interval @f$[a, b]@f$ must contain the root
     * @f$x@f$. The calculations are done with an approximate relative
     * precision `tol`. Returns @f$x@f$ such that @f$f(x) = 0@f$.
     *  @param a            left endpoint of initial interval
     *  @param b            right endpoint of initial interval
     *  @param f            the function which is evaluated
     *  @param tol          accuracy goal
     *  @return the root @f$x@f$
     */
    public static BigDecimal brentDekker (BigDecimal a, BigDecimal b,
                                          MathFunctionBD f, BigDecimal tol) {
        final BigDecimal EPS = new BigDecimal("0.5E-15");
        final int MAXITER = 120;    // Maximum number of iterations
        BigDecimal c, d, e;
        BigDecimal fa, fb, fc;
        final boolean DEBUG = false;

        // Special case I = [b, a]
        if (b.compareTo(a) < 0) {
            BigDecimal ctemp = a;
            a = b;
            b = ctemp;
        }

        // Initialization
        fa = f.evaluate (a);
        if (BigDecimalUtils.compare(fa.abs(), MINVAL, 0) < 0)
            return a;

        fb = f.evaluate (b);
        if (BigDecimalUtils.compare(fb.abs(), MINVAL, 0) <= 0)
            return b;

        c = a;
        fc = fa;
        d = e = b.subtract(a);
        tol = tol.add(EPS).add(DBL_EPSILON); // in case tol is too small

        if (fc.abs().compareTo(fb.abs()) < 0) {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }

        int i;

        for (i = 0; i < MAXITER; i++) {
            BigDecimal s, p, q, r;
            BigDecimal tol1 = tol.add(BigDecimalUtils.FOUR.multiply(DBL_EPSILON.multiply(b.abs())));
            BigDecimal xm = BigDecimalUtils.HALF.multiply(c.subtract(b));
            if (DEBUG) {
                BigDecimal err = fa.subtract(fb).abs();
                System.out.printf("[a, b] = [%g, %g]   fa = %g,   fb = %g   |fa - fb| = %.2g%n",
                        a, b, fa, fb, err);
            }

            if (fb.abs().compareTo(MINVAL) <= 0) {
                return b;
            }
            if (xm.abs().compareTo(tol1) <= 0) {
                if (b.abs().compareTo(MINVAL) > 0)
                    return b;
                else
                    return BigDecimal.ZERO;
            }

            if ((e.abs().compareTo(tol1) >= 0) && (fa.abs().compareTo(fb.abs()) > 0)) {
                if (a.compareTo(c) != 0) {
                    // Inverse quadratic interpolation
                    q = fa.divide(fc, MathContext.DECIMAL64);
                    r = fb.divide(fc, MathContext.DECIMAL64);
                    s = fb.divide(fa, MathContext.DECIMAL64);
                    //s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
                    p = s.multiply( BigDecimalUtils.TWO.multiply(xm).multiply(q).multiply(q.subtract(r))).subtract(b.subtract(a).multiply(r.subtract(BigDecimal.ONE)));
                    //(q - 1.0) * (r - 1.0) * (s - 1.0);
                    q = q.subtract(BigDecimal.ONE).multiply(r.subtract(BigDecimal.ONE)).multiply(s.subtract(BigDecimal.ONE));
                } else {
                    // Linear interpolation
                    s = fb.divide(fa, MathContext.DECIMAL64);
                    p = BigDecimalUtils.TWO.multiply(xm).multiply(s);
                    q = BigDecimal.ONE.subtract(s);
                }

                // Adjust signs
                if (p.compareTo(BigDecimal.ZERO) > 0)
                    q = q.abs();
                p = p.abs();

                // Is interpolation acceptable ?
                //if (((2.0 * p) >= (3.0 * xm * q - Math.abs (tol1 * q)))
                //        || (p >= Math.abs (0.5 * e * q))) {
                if (BigDecimalUtils.TWO.multiply(p).compareTo(BigDecimalUtils.THREE.multiply(xm).multiply(q).subtract(tol1.max(q))) >= 0 ||
                        p.compareTo(BigDecimalUtils.HALF.multiply(e).multiply(q).abs()) >= 0) {

                    d = xm;
                    e = d;
                } else {
                    e = d;
                    d = p.divide(q, MathContext.DECIMAL64);
                }
            } else {
                // Bisection necessary
                d = xm;
                e = d;
            }

            a = b;
            fa = fb;
            if (d.abs().compareTo(tol1) > 0)
                b = b.add(d);
            else if (xm.compareTo(BigDecimal.ZERO) < 0)
                b = b.subtract(tol1);
            else
                b = b.add(tol1);
            fb = f.evaluate (b);
            if (fb.multiply(BigDecimal.valueOf(fc.signum())).compareTo(BigDecimal.ZERO) > 0) {
                c = a;
                fc = fa;
                d = e = b.subtract(a);
            } else {
                a = b;
                b = c;
                c = a;
                fa = fb;
                fb = fc;
                fc = fa;
            }
        }

        if (i >= MAXITER)
            System.err.println(" WARNING:  root finding does not converge");
        return b;
    }

    /**
     * Computes a root @f$x@f$ of the function in `f` using the *bisection*
     * method. The interval @f$[a, b]@f$ must contain the root @f$x@f$. The
     * calculations are done with an approximate relative precision `tol`.
     * Returns @f$x@f$ such that @f$f(x) = 0@f$.
     *  @param a            left endpoint of initial interval
     *  @param b            right endpoint of initial interval
     *  @param f            the function which is evaluated
     *  @param tol          accuracy goal
     *  @return the root @f$x@f$
     */

    public static BigDecimal bisection (BigDecimal a, BigDecimal b,
                                    MathFunctionBD f, BigDecimal tol) {
        // Case I = [b, a]
        if (b.compareTo(a) < 0) {
            BigDecimal ctemp = a;
            a = b;
            b = ctemp;
        }
        BigDecimal xa = a;
        BigDecimal xb = b;
        BigDecimal yb = f.evaluate (b);
        // do preliminary checks on the bounds
        if (yb.abs().compareTo(MINVAL) <= 0)
            return b;
        BigDecimal ya = f.evaluate (a);
        if (ya.abs().compareTo(MINVAL) <= 0)
            return a;

        BigDecimal x = BigDecimal.ZERO, y = BigDecimal.ZERO;
        final int MAXITER = 1200;   // Maximum number of iterations
        final boolean DEBUG = false;
        tol = tol.add(DBL_EPSILON); // in case tol is too small

        if (DEBUG)
            System.out.println
                    ("\niter              xa                   xb              f(x)");

        boolean fini = false;
        int i = 0;
        while (!fini) {
            x = xa.add(xb).divide(BigDecimalUtils.TWO, MathContext.DECIMAL64);
            y = f.evaluate (x);
            if ((y.abs().compareTo(MINVAL) <= 0) ||
                    (xb.subtract(xa).abs().compareTo(tol.multiply(x.abs())) <= 0) ||
                    (xb.subtract(xa).abs()).compareTo(MINVAL) <= 0) {
                if (x.abs().compareTo(MINVAL) > 0)
                    return x;
                else
                    return BigDecimal.ZERO;
            }
            if (y.multiply(ya).compareTo(BigDecimal.ZERO) < 0)
                xb = x;
            else
                xa = x;
            ++i;
            if (DEBUG)
                System.out.printf("%3d    %18.12g     %18.12g    %14.4g%n",
                        i, xa, xb, y);
            if (i > MAXITER) {
                System.out.println ("***** bisection:  SEARCH DOES NOT CONVERGE");
                fini = true;
            }
        }
        return x;
    }

}
