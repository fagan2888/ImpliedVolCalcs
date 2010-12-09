package uk.co.richardbounds;

import org.junit.Test;

import static org.junit.Assert.*;

import uk.co.richardbounds.ImpliedVol.PutCall;

public class ImpliedVolTest {

	@Test
	public void shouldImplyVol() {
		ImpliedVol v = new ImpliedVol();
		double vol = v.gImpliedVolatilityNR(PutCall.PUT, 65, 66, 0.5, 0.1, 0.1, 3, 1e-12); //S, x, T, r, b, cm, epsilon)
		assertEquals(0.2229, vol, 1e-4);
	}
	
}
