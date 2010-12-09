package uk.co.richardbounds;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.NormalDistributionImpl;

public class ImpliedVol {

	enum PutCall {
		PUT, CALL
	};

	/**
	 * Generalised Black Scholes eqn from Haug "Complete guide to option pricing
	 * formulas pg.8".
	 * 
	 * @param putCall
	 * @param s
	 *            stock price
	 * @param x
	 *            strike
	 * @param t
	 *            time to expiry (years)
	 * @param r
	 *            risk-free rate
	 * @param b
	 *            b=r: Black 73 model b=r-q Merton 73 model with dividend q
	 *            b=0 Black 76 futures model 
	 *            b=0 and r=0 Asay 1982 margined futures option model 
	 *            b=r-rf Garman & Kohlhagen 1983 currency option model
	 * @param v
	 *            volatility of the underlying stock price.
	 * @return Option price
	 */
	double gBlackScholes(PutCall putCall, double s, double x, double t,
			double r, double b, double v) {

		double d1, d2;

		d1 = (Math.log(s / x) + (b + Math.pow(v, 2) / 2) * t)
				/ (v * Math.sqrt(t));
		d2 = d1 - v * Math.sqrt(t);

		switch (putCall) {
		case CALL:
			return s * Math.exp((b - r) * t) * cnd(d1) - x
					* Math.exp(-1 * r * t) * cnd(d2);
		case PUT:
			return x * Math.exp(-1 * r * t) * cnd(-1 * d2) - s
					* Math.exp((b - r) * t) * cnd(-1 * d1);
		default:

		}
		return Double.NaN;
	}

	/**
	 *  Vega for the generalized Black and Scholes formula,
	 *  from Haug "Complete guide to option pricing formulas" 
	 * @param s Underlying price.
	 * @param x Strike price.
	 * @param t Time to expiry (years)
	 * @param r Risk free rate.
	 * @param b 
	 *            b=r: Black 73 model b=r-q Merton 73 model with dividend q
	 *            b=0 Black 76 futures model 
	 *            b=0 and r=0 Asay 1982 margined futures option model 
	 *            b=r-rf Garman & Kohlhagen 1983 currency option model.
	 * @param v volatility of the underlying.
	 * @return Option vega.
	 */
	double gVega(double s, double x, double t, double r, double b, double v) {
		double d1 = (Math.log(s / x) + (b + Math.pow(v, 2) / 2) * t)
				/ (v * Math.sqrt(t));
		return s * Math.exp((b - r) * t) * nd(d1) * Math.sqrt(t);
	}

	/**
	 * Newton-Raphson implied vol solver from Haug "Complete guide to option pricing
	 * formulas pg.454".
	 * @param putCall
	 * @param s Stock price.
	 * @param x Strike.
	 * @param t Time to expiry (years)
	 * @param r risk-free rate.
	 * @param b
	 *            b=r: Black 73 model b=r-q Merton 73 model with dividend q
	 *            b=0 Black 76 futures model 
	 *            b=0 and r=0 Asay 1982 margined futures option model 
	 *            b=r-rf Garman & Kohlhagen 1983 currency option model.
	 * @param cm Market price of the option.
	 * @param epsilon
	 * @return
	 */
	double gImpliedVolatilityNR(PutCall putCall, double s, double x,
			double t, double r, double b, double cm, double epsilon) {

		double vi, ci;
		double vegai;
		double minDiff;

		// Manaster and Koehler seed value (vi)
		vi = Math.sqrt(Math.abs(Math.log(s / x) + r * t) * 2 / t);
		ci = gBlackScholes(putCall, s, x, t, r, b, vi);
		vegai = gVega(s, x, t, r, b, vi);
		minDiff = Math.abs(cm - ci);

		while (Math.abs(cm - ci) >= epsilon && Math.abs(cm - ci) <= minDiff) {
			vi = vi - (ci - cm) / vegai;
			ci = gBlackScholes(putCall, s, x, t, r, b, vi);
			vegai = gVega(s, x, t, r, b, vi);
			minDiff = Math.abs(cm - ci);
		}

		if (Math.abs(cm - ci) < epsilon) {
			return vi;
		} else {
			return Double.NaN;
		}

	}

	/**
	 * Normal distribution function.
	 * 
	 * @param x
	 * @return
	 */
	double nd(double x) {
		return new NormalDistributionImpl().density(x);
	}

	/**
	 * Cumulative Normal distribution function.
	 * 
	 * @param x
	 * @return
	 */
	double cnd(double x) {
		try {
			return new NormalDistributionImpl().cumulativeProbability(x);
		} catch (MathException e) {
			throw new RuntimeException(e);
		}
	}

}
