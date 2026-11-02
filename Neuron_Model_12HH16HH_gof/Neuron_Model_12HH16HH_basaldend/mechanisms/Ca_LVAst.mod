:Comment : LVA ca channel. Note: mtau is an approximation from the plots
:Reference : :		Avery and Johnston 1996, tau from Randall 1997
:Comment: shifted by -10 mv to correct for junction potential
:Comment: corrected rates using q10 = 2.3, target temperature 34, orginal 21
:
: TF 2026-08-21 K1422E: added actshift / inactshift RANGE parameters.
:   BOTH DEFAULT TO 0, so with no explicit assignment this file is numerically
:   identical to the original (backed up as Ca_LVAst.mod.pre_shift_backup).
:
:   As published here the curves sit at
:       activation   V1/2 = -40.0 mV (slope 6)
:       inactivation V1/2 = -90.0 mV (slope 6.4)
:   which leaves the channel ~94% inactivated at the model's -72 mV resting
:   potential, so it carries essentially no window current and a T-type block
:   is unmeasurable. The Hallermann/Ben-Shalom K1422E model instead used the
:   Destexhe/Huguenard `it2` current shifted to
:       activation   V1/2 = -55 mV (slope 7.4)
:       inactivation V1/2 = -75 mV (slope 5.0)
:   i.e. actshift = -15, inactshift = +15.
:
:   V1/2_act   = -40 + actshift
:   V1/2_inact = -90 + inactshift
:   Time constants shift with their own curve.

NEURON	{
	SUFFIX Ca_LVAst
	USEION ca READ eca WRITE ica
	NONSPECIFIC_CURRENT i_nc
	RANGE gCa_LVAstbar, gCa_LVAst, ica, i_nc, itot, pool_frac
	RANGE actshift, inactshift, mslope, hslope
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gCa_LVAstbar = 0.00001 (S/cm2)
	actshift   = 0 (mV)  : moves activation   V1/2 to -40 + actshift
	inactshift = 0 (mV)  : moves inactivation V1/2 to -90 + inactshift
: ---- pool_frac: nanodomain separation of T-type Ca2+ from SK ---------------
: CaDynamics_E2 lumps EVERY Ca2+ source in a segment into ONE `cai`, which then
: drives SK_E2. So Ca_LVAst is wired directly into SK even though SK in cortical
: pyramidal cells is coupled in nanodomains predominantly to HVA / L-type
: channels. The artefact is large: blocking T-type unloads cai, shuts SK off and
: INCREASES firing -- the opposite of what ethosuximide does clinically.
:
: pool_frac is the fraction of the T-type current that reaches the shared cai
: pool. The TOTAL membrane current is unchanged (ica + i_nc = itot at every
: timestep) -- only its visibility to CaDynamics_E2 -> SK_E2 changes.
:   pool_frac = 1 (default)  identical to the published model
:   pool_frac = 0            T-type still carries its full current but is
:                            invisible to the SK-driving Ca2+ pool
	pool_frac  = 1       : fraction of I_CaT that feeds the shared cai pool
	mslope     = 6   (mV): activation slope   (published 6)
	hslope     = 6.4 (mV): inactivation slope (published 6.4)
	: Old Hallermann/Ben-Shalom `it2` equals
	:     actshift = -15, mslope = 7.4, inactshift = +15, hslope = 5.0
}

ASSIGNED	{
	v	(mV)
	eca	(mV)
	ica	(mA/cm2)
	gCa_LVAst	(S/cm2)
	i_nc	(mA/cm2)
	itot	(mA/cm2)
	mInf
	mTau
	hInf
	hTau
}

STATE	{
	m
	h
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	gCa_LVAst = gCa_LVAstbar*m*m*h
	itot = gCa_LVAst*(v-eca)
	ica  = itot*pool_frac
	i_nc = itot*(1-pool_frac)
}

DERIVATIVE states	{
	rates()
	m' = (mInf-m)/mTau
	h' = (hInf-h)/hTau
}

INITIAL{
	rates()
	m = mInf
	h = hInf
}

PROCEDURE rates(){
  LOCAL qt
  qt = 2.3^((34-21)/10)

	UNITSOFF
		: V1/2 act = -40 + actshift ; V1/2 inact = -90 + inactshift
		mInf = 1.0000/(1+ exp((v - (-40.000 + actshift))/(-mslope)))
		mTau = (5.0000 + 20.0000/(1+exp((v - (-35.000 + actshift))/5)))/qt
		hInf = 1.0000/(1+ exp((v - (-90.000 + inactshift))/hslope))
		hTau = (20.0000 + 50.0000/(1+exp((v - (-50.000 + inactshift))/7)))/qt
	UNITSON
}
