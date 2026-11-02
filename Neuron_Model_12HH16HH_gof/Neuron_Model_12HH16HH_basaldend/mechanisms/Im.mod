:Reference : :		Adams et al. 1982 - M-currents and other potassium currents in bullfrog sympathetic neurones
:Comment: corrected rates using q10 = 2.3, target temperature 34, orginal 21

NEURON	{
	SUFFIX Im
	USEION k READ ek WRITE ik
	RANGE gImbar, gIm, ik, vshift, taufac
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gImbar = 0.00001 (S/cm2) 
	: vshift moves the activation curve depolarised. 0 = published Adams
	: kinetics, i.e. this file is a strict no-op unless model_config sets it.
	vshift = 0 (mV)
	: taufac DIVIDES mTau. 1 = published Adams kinetics (no-op).
	: WHY IT EXISTS. The Scn2a KO's firing advantage in this model comes from a
	: MISSING Im brake, and Im's mTau here is ~150 ms at threshold voltages, so
	: the brake only accumulates across a spike TRAIN. That is why WT and KO
	: overlap exactly up to 175 pA (0-2 spikes, no train, no brake difference)
	: and separate only once repetitive firing begins. Spratt Fig 1B shows the
	: KO above WT across the measured range. Speeding Im up moves the divergence
	: point to lower currents without changing the brake's STRENGTH (mInf is
	: untouched), so gbar and vshift keep their meaning.
	taufac = 1 (1)
}

ASSIGNED	{
	v	(mV)
	ek	(mV)
	ik	(mA/cm2)
	gIm	(S/cm2)
	mInf
	mTau
	mAlpha
	mBeta
}

STATE	{ 
	m
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	gIm = gImbar*m
	ik = gIm*(v-ek)
}

DERIVATIVE states	{
	rates()
	m' = (mInf-m)/mTau
}

INITIAL{
	rates()
	m = mInf
}

PROCEDURE rates(){
  LOCAL qt
  qt = 2.3^((34-21)/10)

	UNITSOFF
		mAlpha = 3.3e-3*exp(2.5*0.04*(v - -35 - vshift))
		mBeta = 3.3e-3*exp(-2.5*0.04*(v - -35 - vshift))
		mInf = mAlpha/(mAlpha + mBeta)
		mTau = (1/(mAlpha + mBeta))/qt/taufac
	UNITSON
}
