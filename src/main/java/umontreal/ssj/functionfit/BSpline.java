/*
 * Class:        BSpline
 * Description:  B-spline
 * Environment:  Java
 * Software:     SSJ 
 * Copyright (C) 2001  Pierre L'Ecuyer and Universite de Montreal
 * Organization: DIRO, Universite de Montreal
 * @author       
 * @since
 *
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
package umontreal.ssj.functionfit;

import umontreal.ssj.functions.*;
import umontreal.ssj.util.BigDecimalUtils;
import umontreal.ssj.util.Misc;
import umontreal.ssj.util.RootFinder;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.DoubleFactory2D;
import umontreal.ssj.util.RootFinderBD;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.Arrays;
import java.util.Random;
import java.io.*;
import java.util.zip.ZipError;

/**
 * Represents a B-spline with control points at @f$(X_i, Y_i)@f$. Let
 * @f$\mathbf{P_i}=(X_i, Y_i)@f$, for @f$i=0,…,n-1@f$, be a *control point*
 * and let @f$t_j@f$, for @f$j=0,…,m-1@f$ be a *knot*. A B-spline
 * @cite mDEB78a&thinsp; of degree @f$p=m-n-1@f$ is a parametric curve
 * defined as
 * @f[
 *   \mathbf{P(t)} = \sum_{i=0}^{n-1} N_{i, p}(t) \mathbf{P_i},\mbox{ for }t_p\le t\le t_{m-p-1}.
 * @f]
 * Here,
 * @f{align*}{
 *    N_{i, p}(t) 
 *    & 
 *   =
 *    \frac{t-t_i}{t_{i+p} - t_i}N_{i, p-1}(t) + \frac{t_{i+p+1} - t}{t_{i+p+1} - t_{i+1}}N_{i+1, p-1}(t)
 *    \\ 
 *    N_{i, 0}(t) 
 *    & 
 *   =
 *    \left\{
 *   \begin{array}{ll}
 *    1
 *    & 
 *   \mbox{ for }t_i\le t\le t_{i+1},
 *    \\ 
 *   0
 *   \mbox{ elsewhere}. 
 *   \end{array}
 *   \right.
 * @f}
 * This class provides methods to evaluate @f$\mathbf{P(t)}=(X(t), Y(t))@f$
 * at any value of @f$t@f$, for a B-spline of any degree @f$p\ge1@f$. Note
 * that the  #evaluate(double) method of this class can be slow, since it
 * uses a root finder to determine the value of @f$t^*@f$ for which
 * @f$X(t^*)=x@f$ before it computes @f$Y(t^*)@f$.
 *
 * <div class="SSJ-bigskip"></div>
 */
public class BSpline implements MathFunctionBD

,
MathFunctionWithIntegral, MathFunctionWithDerivative, MathFunctionWithFirstDerivative{
   // Formula taken from http://www.ibiblio.org/e-notes/Splines/Basis.htm
   // and http://en.wikipedia.org/wiki/De_Boor_algorithm
   private BigDecimal[] x;     //x original
   private BigDecimal[] y;     //y original

   private int degree;

   //variables sur lesquelles on travaille
   private BigDecimal[] myX;
   private BigDecimal[] myY;
   private BigDecimal[] knots;

   private static BigDecimal[][] initBDMatrix (int i, int j) {
      BigDecimal[][] ret = new BigDecimal[i][j];
      for (BigDecimal[] bigDecimals : ret) Arrays.fill(bigDecimals, BigDecimal.ZERO);

      return ret;
   }

   private static BigDecimal[] emptyArray(int length) {
      BigDecimal[] ret = new BigDecimal[length];
      Arrays.fill(ret, BigDecimal.ZERO);
      return ret;
   }

   /**
    * Constructs a new uniform B-spline of degree `degree` with control
    * points at (<tt>x[i]</tt>, <tt>y[i]</tt>). The knots of the resulting
    * B-spline are set uniformly from `x[0]` to `x[n-1]`.
    *  @param x            the values of @f$X@f$.
    *  @param y            the values of @f$Y@f$.
    *  @param degree       the degree of the B-spline.
    */
   public BSpline (final BigDecimal[] x, final BigDecimal[] y, final int degree) {
      if (x.length != y.length)
         throw new IllegalArgumentException("The arrays x and y must share the same length");
      this.degree = degree;
      this.x = x.clone();
      this.y = y.clone();
      init(x, y, null);
   }

   /**
    * Constructs a new uniform B-spline with control points at
    * (<tt>x[i]</tt>, <tt>y[i]</tt>), and knot vector given by the array
    * `knots`.
    *  @param x            the values of @f$X@f$.
    *  @param y            the values of @f$Y@f$.
    *  @param knots        the knots of the B-spline.
    */
   public BSpline (final BigDecimal[] x, final BigDecimal[] y, final BigDecimal[] knots) {
      if (x.length != y.length)
         throw new IllegalArgumentException("The arrays x and y must share the same length");
      if (knots.length < x.length + 1)
         throw new IllegalArgumentException("The number of knots must be at least n+1");

      this.x = x.clone();
      this.y = y.clone();
      this.knots = knots.clone();
      init(x, y, knots);
   }

   /**
    * Returns the @f$X_i@f$ coordinates for this spline.
    *  @return the @f$X_i@f$ coordinates.
    */
   public BigDecimal[] getX() {
      return myX.clone ();
   }

   /**
    * Returns the @f$Y_i@f$ coordinates for this spline.
    *  @return the @f$Y_i@f$ coordinates.
    */
   public BigDecimal[] getY() {
      return myY.clone ();
   }

   /**
    * Returns the knot maximal value.
    *  @return the @f$Y_i@f$ coordinates.
    */
   public BigDecimal getMaxKnot() {
      return knots[knots.length-1];
   }

   /**
    * Returns the knot minimal value.
    *  @return the @f$Y_i@f$ coordinates.
    */
   public BigDecimal getMinKnot() {
      return knots[0];
   }

   /**
    * Returns an array containing the knot vector @f$(t_0, t_{m-1})@f$.
    *  @return the knot vector.
    */
   public BigDecimal[] getKnots() {
      return knots.clone ();
   }

   /**
    * Returns a B-spline curve of degree `degree` interpolating the
    * @f$(x_i, y_i)@f$ points @cite mDEB78a&thinsp;. This method uses the
    * uniformly spaced method for interpolating points with a B-spline
    * curve, and a uniformed clamped knot vector, as described in
    * [http://www.cs.mtu.edu/~shene/COURSES/cs3621/NOTES/](http://www.cs.mtu.edu/~shene/COURSES/cs3621/NOTES/).
    *  @param x            the values of @f$X@f$.
    *  @param y            the values of @f$Y@f$.
    *  @param degree       the degree of the B-spline.
    *  @return the B-spline curve.
    */
   public static BSpline createInterpBSpline (BigDecimal[] x, BigDecimal[] y,
                                              int degree) {
      if (x.length != y.length)
         throw new IllegalArgumentException("The arrays x and y must share the same length");
      if (x.length <= degree)
         throw new IllegalArgumentException("The arrays length must be greater than degree");

      int n = x.length-1;
      //compute t : parameters vector uniformly from 0 to 1
      BigDecimal[] t = emptyArray(x.length);
      for(int i = 0; i<t.length; i++) {
         t[i] = BigDecimal.valueOf((double)i/(t.length-1));
      }

      //compute U : clamped knots vector uniformly from 0 to 1
      BigDecimal[] U = emptyArray(x.length + degree + 1);
      int m = U.length-1;
      for(int i =0; i<=degree; i++)
         U[i] = BigDecimal.ZERO;
      for(int i =1; i<x.length-degree; i++)
         U[i+degree] = BigDecimal.valueOf((double)i/(x.length-degree));
      for(int i = U.length-1-degree; i<U.length; i++)
         U[i] = BigDecimal.ONE;


      //compute matrix N : made of BSpline coefficients
      BigDecimal [][] N = initBDMatrix(x.length, x.length);
      for(int i = 0; i<x.length; i++) {
            N[i] = computeN(U, degree, t[i], x.length);
      }

      //initialize D : initial points matrix
      BigDecimal [][] D = initBDMatrix(x.length, 2);
      for(int i =0; i<x.length; i++) {
         D[i][0] = x[i];
         D[i][1] = y[i];
      }

      //solve the linear equation system using colt library
      DoubleMatrix2D coltN = DoubleFactory2D.dense.make(BigDecimalUtils.toDoubleMatrix(N));
      DoubleMatrix2D coltD = DoubleFactory2D.dense.make(BigDecimalUtils.toDoubleMatrix(D));
      DoubleMatrix2D coltP;
      coltP = Algebra.ZERO.solve(coltN, coltD);

      return new BSpline(BigDecimalUtils.toDecimalArray(coltP.viewColumn(0).toArray(), MathContext.DECIMAL64),
              BigDecimalUtils.toDecimalArray(coltP.viewColumn(1).toArray(), MathContext.DECIMAL64), U);
   }
   
   
   /**
    * Returns a B-spline curve of degree `degree` smoothing @f$(x_i,
    * y_i)@f$, for @f$i=0,…,n@f$ points. The precision depends on the
    * parameter @f$hp1@f$: @f$1 \le\mathtt{degree} \le hp1<n@f$, which
    * represents the number of control points used by the new B-spline
    * curve, minimizing the quadratic error
    * @f[
    *   L = \sum_{i=0}^n\left( \frac{Y_i - S_i(X_i)}{W_i}\right)^2.
    * @f]
    * This method uses the uniformly spaced method for interpolating
    * points with a B-spline curve and a uniformed clamped knot vector, as
    * described in
    * [http://www.cs.mtu.edu/~shene/COURSES/cs3621/NOTES/](http://www.cs.mtu.edu/~shene/COURSES/cs3621/NOTES/).
    *  @param x            the values of @f$X@f$.
    *  @param y            the values of @f$Y@f$.
    *  @param degree       the degree of the B-spline.
    *  @param hp1          the desired number of control points.
    *  @return the B-spline curve.
    */
   public static BSpline createApproxBSpline (BigDecimal[] x, BigDecimal[] y,
                                              int degree, int hp1) {
      if (x.length != y.length)
         throw new IllegalArgumentException("The arrays x and y must share the same length");
      if (x.length <= degree)
         throw new IllegalArgumentException("The arrays length must be greater than degree");

      int h = hp1 - 1;
      int n = x.length-1;
      
      //compute t : parameters vector uniformly from 0 to 1
      BigDecimal[] t = emptyArray(x.length);
      for(int i = 0; i < t.length; i++) {
         t[i] = BigDecimal.valueOf((double)i / n);
      }

      //compute U : clamped knots vector uniformly from 0 to 1
      int m = h + degree + 1;
      BigDecimal U[] = emptyArray(m + 1);
      for(int i = 0; i <= degree; i++)
         U[i] = BigDecimal.ZERO;
      for(int i = 1; i < hp1 - degree; i++)
         U[i+degree] = BigDecimal.valueOf((double)i/(hp1 - degree));
      for(int i = m-degree; i<= m; i++)
         U[i] = BigDecimal.ONE;

      
      //compute matrix N : composed of BSpline coefficients
      BigDecimal [][] N = initBDMatrix(n+1, h+1);
      for(int i = 0; i < N.length; i++) {
         N[i] = computeN(U, degree, t[i], h+1);
      }

      //initialize D : initial points matrix
      BigDecimal [][] D = initBDMatrix(x.length, 2);
      for(int i = 0; i < x.length; i++) {
         D[i][0] = x[i];
         D[i][1] = y[i];
      }

      //compute Q :
      BigDecimal[][] tempQ = initBDMatrix(x.length, 2);
      for(int k = 1; k < n; k++) {
         tempQ[k][0] = D[k][0].subtract(N[k][0].multiply(D[0][0])).subtract(N[k][h].multiply(D[D.length-1][0]));
         tempQ[k][1] = D[k][1].subtract(N[k][0].multiply(D[0][1])).subtract(N[k][h].multiply(D[D.length-1][1]));
      }
      BigDecimal[][] Q = initBDMatrix(h-1, 2);
      for(int i = 1; i < h; i++) {
         for(int k = 1; k < n; k++) {
            Q[i-1][0] = Q[i-1][0].add(N[k][i].multiply(tempQ[k][0]));
            Q[i-1][1] = Q[i-1][1].add(N[k][i].multiply(tempQ[k][1]));
         }
      }
      
      // compute N matrix for computation:
      BigDecimal [][] N2 = initBDMatrix(n-1, h-1);
      for(int i = 0; i < N2.length; i++) {
         if (h - 1 >= 0) System.arraycopy(N[i + 1], 1, N2[i], 0, h - 1);
      }

      //solve the linear equation system using colt library
      DoubleMatrix2D coltQ = DoubleFactory2D.dense.make(BigDecimalUtils.toDoubleMatrix(Q));
      DoubleMatrix2D coltN = DoubleFactory2D.dense.make(BigDecimalUtils.toDoubleMatrix(N2));
      DoubleMatrix2D coltM = Algebra.ZERO.mult(Algebra.ZERO.transpose(coltN), coltN);
      DoubleMatrix2D coltP = Algebra.ZERO.solve(coltM, coltQ);
      BigDecimal[] pxTemp = BigDecimalUtils.toDecimalArray(coltP.viewColumn(0).toArray(), MathContext.DECIMAL64);
      BigDecimal[] pyTemp = BigDecimalUtils.toDecimalArray(coltP.viewColumn(1).toArray(), MathContext.DECIMAL64);
      BigDecimal[] px = emptyArray(hp1);
      BigDecimal[] py = emptyArray(hp1);
      px[0] = D[0][0];
      py[0] = D[0][1];
      px[h] = D[D.length-1][0];
      py[h] = D[D.length-1][1];
      for(int i = 0; i < pxTemp.length; i++) {
         px[i+1] = pxTemp[i];
         py[i+1] = pyTemp[i];
      }

      return new BSpline(px, py, U);
      // return new BSpline(px, py, degree);
   }
   
   
   /**
    * Returns the derivative B-spline object of the current variable.
    * Using this function and the returned object, instead of the
    * `derivative` method, is strongly recommended if one wants to compute
    * many derivative values.
    *  @return the derivative B-spline of the current variable.
    */
   public BSpline derivativeBSpline() {
      BigDecimal[] xTemp = emptyArray(this.myX.length - 1);
      BigDecimal[] yTemp = emptyArray(this.myY.length - 1);
      BigDecimal deg = BigDecimal.valueOf(degree);
      for(int i = 0; i<xTemp.length; i++) {
         xTemp[i] = myX[i+1].subtract(myX[i]).multiply(deg).divide(knots[i+degree+1].subtract(knots[i+1]), MathContext.DECIMAL64);
         yTemp[i] = myY[i+1].subtract(myY[i]).multiply(deg).divide(knots[i+degree+1].subtract(knots[i+1]), MathContext.DECIMAL64);
      }

      BigDecimal [] newKnots = emptyArray(knots.length - 2);
      System.arraycopy(knots, 1, newKnots, 0, newKnots.length);

      //tri pas optimise du tout
      BigDecimal[] xTemp2 = emptyArray(this.myX.length-1);
      BigDecimal[] yTemp2 = emptyArray(this.myY.length-1);
      for(int i = 0; i<xTemp.length; i++) {
         int k=0;
         for (BigDecimal bigDecimal : xTemp) {
            if (xTemp[i].compareTo(bigDecimal) > 0)
               k++;
         }
         while(!xTemp2[k].equals(BigDecimal.ZERO))
            k++;
         xTemp2[k] = xTemp[i];
         yTemp2[k] = yTemp[i];
      }

      return new BSpline(xTemp2, yTemp2, newKnots);
   }

   /**
    * Returns the @f$i@f$th derivative B-spline object of the current
    * variable; @f$i@f$ must be less than the degree of the original
    * B-spline. Using this function and the returned object, instead of
    * the `derivative` method, is strongly recommended if one wants to
    * compute many derivative values.
    *  @param i            the degree of the derivative.
    *  @return the ith derivative.
    */
   public BSpline derivativeBSpline (int i) {
      BSpline bs = this;
      while(i > 0) {
         i--;
         bs = bs.derivativeBSpline();
      }
      return bs;
   }

   public double evaluate(final double u) {
      return evaluate(BigDecimal.valueOf(u)).doubleValue();
   }

   public BigDecimal evaluate(final BigDecimal u) {
      final MathFunctionBD xFunction = new MathFunctionBD () {
         @Override
         public BigDecimal evaluate(BigDecimal x) {
            return evalX(x).subtract(u);
         }
      };
      // brentDekker may be unstable; using bisection instead
      // final double t = RootFinder.brentDekker (0, 1, xFunction, 1e-6);
      final BigDecimal t = RootFinderBD.bisection(BigDecimal.ZERO, BigDecimal.ONE, xFunction, new BigDecimal("1e-6"));
      return evalY(t);
   }

   public double integral (double a, double b) {
      return MathFunctionUtil.simpsonIntegral (this, a, b, 500);
   }

   public double derivative(double u) {
      return derivativeBSpline().evaluate(u);
   }

   public double derivative(double u, int n) {
      return derivativeBSpline(n).evaluate(u);
   }

   private void init(BigDecimal[] x, BigDecimal[] y, BigDecimal [] initialKnots) {
      if(initialKnots == null) {
         //Cree un vecteur de noeud uniforme entre 0 et 1
         knots = emptyArray(x.length + degree + 1);

         for(int i = degree; i < this.knots.length-degree; i++)
            this.knots[i]= BigDecimal.valueOf((double)(i-degree)/(knots.length - (2.0*degree) -1));

         if (this.knots.length - (this.knots.length - degree) >= 0)
            System.arraycopy(this.knots, this.knots.length - degree - 1, this.knots, this.knots.length - degree, this.knots.length - (this.knots.length - degree));

         if (degree >= 0) System.arraycopy(this.knots, 1, this.knots, 0, degree);

         // cree notre vecteur interne de Points de controle
         // ici, aucune modification a faire sur les tableaux originaux
         myX = x;
         myY = y;
      }
      else {
         degree = initialKnots.length - x.length -1;

      // on adapte le tableau des noeuds a notre algorithme
      // le tableau knot necessite d'avoir degree fois la meme valeur en debut et en fin de tableau
      // on adapte la taille des tableau X et Y en consequence afin de continuer a respecter la condition :
      // x.length + degree + 1 = this.knots.length
      // Cette modification n'influence pas le resultat et permet de fermer notre courbe

         //Compute the number of values that need to be added
         int iBorneInf = 1, iBorneSup = initialKnots.length-2;
         // while(initialKnots[iBorneInf] == initialKnots[0])

         while (BigDecimalUtils.equals(initialKnots[iBorneInf], initialKnots[0], 10))
            iBorneInf++;
         if (iBorneInf <= degree)
            iBorneInf = degree-iBorneInf+1;
         else
            iBorneInf=0;//on a alors iBorneInf valeurs a rajouter en debut de tableau

         // while(initialKnots[iBorneSup] == initialKnots[initialKnots.length-1])
         while (BigDecimalUtils.equals(initialKnots[iBorneSup], initialKnots[initialKnots.length-1], 10))
            iBorneSup--;
         if (iBorneSup >= initialKnots.length-1-degree)
            iBorneSup = degree+1-(initialKnots.length-1-iBorneSup);
         else
            iBorneSup = 0; //on a alors iBorneSup valeurs a rajouter en fin de tableau

         //add computed values
         this.knots = emptyArray(initialKnots.length + iBorneInf + iBorneSup);
         myX = emptyArray(x.length + iBorneInf + iBorneSup);
         myY = emptyArray(y.length + iBorneInf + iBorneSup);
         for (int i = 0; i < iBorneInf; i++) {
            this.knots[i] = initialKnots[0];
            myX[i] = x[0];
            myY[i] = y[0];
         }
         System.arraycopy(initialKnots, 0, this.knots, iBorneInf, initialKnots.length);
         for(int i = 0; i<x.length; i++) {
            myX[iBorneInf + i] = x[i];
            myY[iBorneInf + i] = y[i];
         }
         for(int i = 0; i<iBorneSup; i++) {
            this.knots[this.knots.length-1 - i] = initialKnots[initialKnots.length-1];
            myX[myX.length-1-i] = x[x.length-1];
            myY[myY.length-1-i] = y[y.length-1];
         }
      }
   }

   public BigDecimal evalX(BigDecimal u) {
      int k = Misc.getTimeInterval (knots, 0, knots.length - 1, u);
      if (k >= myX.length) // set to last control point
         k = myX.length-1;
      
      BigDecimal[][] X = initBDMatrix(degree + 1, myX.length);

      if (k + 1 - (k - degree) >= 0) System.arraycopy(myX, k - degree, X[0], k - degree, k + 1 - (k - degree));

      for(int j=1; j<= degree; j++) {
         for(int i = k-degree+j; i <= k; i++) {
            BigDecimal aij = u.subtract(this.knots[i]).divide(this.knots[i+1+degree-j].subtract(this.knots[i]), MathContext.DECIMAL64);
            X[j][i] = BigDecimal.ONE.subtract(aij).multiply(X[j-1][i-1]).add(aij.multiply(X[j-1][i]));
         }
      }
      return X[degree][k];
   }

   public BigDecimal evalY(BigDecimal u) {
      int k = Misc.getTimeInterval (knots, 0, knots.length - 1, u);
      if (k >= myY.length) // set to last control point
         k = myY.length-1;
      
      BigDecimal[][] Y = initBDMatrix(degree + 1, myX.length);

      if (k + 1 - (k - degree) >= 0) System.arraycopy(myY, k - degree, Y[0], k - degree, k + 1 - (k - degree));

      for(int j=1; j<= degree; j++) {
         for(int i = k-degree+j; i <= k; i++) {
            BigDecimal aij = u.subtract(this.knots[i]).divide(this.knots[i+1+degree-j].subtract(this.knots[i]), MathContext.DECIMAL64);
            Y[j][i] = BigDecimal.ONE.subtract(aij).multiply(Y[j-1][i-1]).add(aij.multiply(Y[j-1][i]));
         }
      }
      return Y[degree][k];
   }

   /**
    * Checks if two doubles {@code a} and {@code b} are equal with tolerance level.
    * 
    * @param a
    * @param b
    * @param tol absolute tolerance level
    * @return 
    */
   private static boolean areEqual(double a, double b, double tol) {
      return Math.abs(a - b) < tol;
   }
   
   private static BigDecimal[] computeN(BigDecimal[] U, int degree, BigDecimal u, int np1) {
      BigDecimal[] N = emptyArray(np1);

      // special cases at bounds
      if (BigDecimalUtils.equals(u, U[0], 10)) {
         N[0] = BigDecimal.ONE;
         return N;
      }
      else if (BigDecimalUtils.equals(u, U[U.length-1], 10)) {
         N[N.length-1] = BigDecimal.ONE;
         return N;
      }

      // find the knot index k such that U[k]<= u < U[k+1]
      int k = Misc.getTimeInterval (U, 0, U.length - 1, u);

      N[k] = BigDecimal.ONE;
      for(int d = 1; d <= degree; d++) {
         N[k-d] = N[k-d+1].multiply(U[k+1].subtract(u)).divide(U[k+1].subtract(U[k-d+1]), MathContext.DECIMAL64);

         for(int i = k-d+1; i<= k-1; i++)
            N[i] = u.subtract(U[i]).divide(U[i+d].subtract(U[i]).multiply(N[i]), MathContext.DECIMAL64).add(
                    U[i+d+1].subtract(u).divide(U[i+d+1].subtract(U[i+1]), MathContext.DECIMAL64).multiply(N[i+1]));

         N[k] = u.subtract(U[k]).divide(U[k+d].subtract(U[k]), MathContext.DECIMAL64).multiply(N[k]);
      }
      return N;
   }

}