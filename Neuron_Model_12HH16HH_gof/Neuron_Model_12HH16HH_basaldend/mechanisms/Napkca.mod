TITLE na12 K1422E -- non-selective cation channel (Na+ / K+ / Ca2+)
: ============================================================================
: SCN2A / Nav1.2 K1422E "generic cation channel"
:
: K1422E sits in the DIII pore loop and converts the DEKA selectivity filter to
: DEEA. The channel still GATES like Nav1.2, but it no longer selects for Na+ --
: it conducts Na+, K+ and Ca2+, so its reversal potential collapses from ~+55 mV
: toward ~0 mV.
:
: PROVENANCE
:   Gating (STATE m,h,s; PROCEDURE trates; FUNCTIONs alpv/alps/bets/trap0 and
:   every rate PARAMETER) is copied VERBATIM from na12HHmut.mod. Nothing about
:   activation, inactivation or slow inactivation has been touched. The ONLY
:   change is how the open conductance is converted into current.
:
:   Replaces Napkca.mod from the Ben-Shalom K1422E repo, which did not compile
:   ("USEION k READ kna WRITE kna" -- kna is not a valid k-ion variable) and
:   which hard-coded an equal 1/3 : 1/3 : 1/3 split.
:
: HOW PERMEATION WORKS
:   na12mut computes one current:      ina = thegna * (v - Ena)
:   napkca  splits the SAME thegna into three parallel ohmic arms, each with its
:   own driving force:
:       ina = thegna * (pna/ptot) * (v - Ena)
:       ik  = thegna * (pk /ptot) * (v - Ek )
:       ica = thegna * (pca/ptot) * (v - Eca)
:   where ptot = pna + pk + pca.
:
:   Total conductance is therefore IDENTICAL to na12mut at the same gbar -- only
:   the ion it is shared among changes. The net reversal potential is the
:   permeability-weighted mean of the three reversals:
:
:       erev = (pna*Ena + pk*Ek + pca*Eca) / (pna + pk + pca)
:
:   Defaults pna:pk:pca = 1 : 0.8 : 0.1 with Ena=55, Ek=-85, Eca=130 give
:       erev = (55 - 68 + 13)/1.9 = 0.0 mV      <-- non-selective cation channel
:   Set pk = pca = 0 to recover exactly the WT-selectivity na12mut behaviour.
:
: REVERSAL POTENTIALS
:   By default the driving forces use the PARAMETERs Ena/Ek/Eca, matching the
:   convention already used by na12/na12mut/na16/na16mut (all of which compute
:   from PARAMETER Ena = 55 and ignore the section's ion ena). This keeps erev
:   exactly predictable and tunable, and is immune to the fact that basal
:   dendrites in biophysics.hoc never get an explicit ek.
:   Set use_ion_erev = 1 to instead read the section's ena/ek/eca, in which case
:   the Ca2+ arm tracks cai through CaDynamics_E2.
:
:   Note: ica is written to the ca ion either way, so it feeds CaDynamics_E2 ->
:   cai -> SK_E2 exactly as a real Ca2+ influx would. Set pca = 0 to decouple.
: ============================================================================

NEURON {
	SUFFIX napkca
	USEION na READ ena WRITE ina
	USEION k  READ ek  WRITE ik
	USEION ca READ eca WRITE ica
	RANGE  gbar, ar2, thegna, ina_ina
	RANGE  pna, pk, pca, Ena, Ek, Eca, use_ion_erev
	RANGE  gna, gk, gca, erev, ik_ik, ica_ica, itot
	RANGE  pers_frac, pers_vhalf, pers_slope, pers_h_vhalf, pers_h_slope
	RANGE  thegna_pers, thegna_trans
	RANGE vhalfs,sh,tha,qa,Ra,Rb,thi1,thi2,qd,qg,mmin,hmin,q10,Rg,qq,Rd,tq,thinf,qinf,vhalfs,a0s,zetas,gms,smax,vvh,vvs
}

PARAMETER {
	: ---- selectivity / permeation (THE K1422E PART) ------------------------
	pna  = 1.0					: relative permeability to Na+ (reference = 1)
	pk   = 0.8					: relative permeability to K+
	pca  = 0.1					: relative permeability to Ca2+
	Ena = 55	(mV)			: Na+ reversal   (matches na12/na12mut/na16)
	Ek  = -85	(mV)			: K+  reversal   (matches biophysics.hoc ek)
	Eca = 130	(mV)			: Ca2+ reversal  (~Nernst at cai=1e-4, cao=2)
	use_ion_erev = 0			: 0 = use Ena/Ek/Eca above; 1 = use ion ena/ek/eca

	: ---- PERSISTENT arm (TF 2026-08-21) -----------------------------------
	: pers_frac DEFAULTS TO 0, so this file is numerically identical to the
	: original (Napkca.mod.pre_pers_backup) unless it is set explicitly.
	:
	: WHY IT EXISTS. thegna = gbar*m^3*h*s is fully TRANSIENT: at the model's
	: -72 mV resting potential m^3*h*s = 1e-12, so the K1422E pore is shut
	: except during the ~0.5 ms AP upstroke, and its Ca2+ arm can therefore
	: never produce a TONIC Ca2+ load. The Hallermann/Ben-Shalom K1422E model
	: implemented the mutant channel as Nap2 -- a PERSISTENT (Magistretti &
	: Alonso 1999) Na channel with 1/3 of its conductance carrying Ca2+ --
	: whose open fraction at -72 mV is 1.3e-2, i.e. 1e10 x larger.
	: This arm reproduces that carrier as an add-on sized by pers_frac.
	pers_frac    = 0				: fraction of gbar in the persistent arm
	pers_vhalf   = -52.6	(mV)	: Magistretti & Alonso 1999 (old Nap2)
	pers_slope   = 4.6		(mV)
	pers_h_vhalf = -48.8	(mV)	: slow inactivation of the persistent arm
	pers_h_slope = 10		(mV)

	: ---- gating: VERBATIM from na12HHmut.mod -------------------------------
	sh   = 8.886047186457889	(mV)
	gbar = 0.01    				(mho/cm2)

	tha  =  -24.155451306086988	(mV)		: v 1/2 for act
	qa   = 5.41					(mV)		: act slope (4.5)
	Ra   = 0.3380714915775742 	(/ms)		: open (v)
	Rb   = 0.09013136340161398 	(/ms)		: close (v)

	thi1  = -60.488477521934875	(mV)		: v 1/2 for inact
	thi2  = -77.41692349310195 	(mV)		: v 1/2 for inact
	qd   = 0.8058343822410788	(mV)	        : inact tau slope
	qg   = 0.6693522946835427    (mV)
	mmin = 0.013671131800210966
	hmin = 0.008420778920829085
	q10 = 2.289601426305275
	Rg   = 0.01854277725353276 	(/ms)		: inact recov (v)
	Rd   = 0.025712394696815438 (/ms)		: inact (v)
	qq   = 10        			(mV)
	tq   = -55      			(mV)

	thinf  = -40.114984963535186  	(mV)		: inact inf slope
	qinf  = 5.760329120353593		(mV)		: inact inf slope

        vhalfs = -33.73363659219147	(mV)		: slow inact.
        a0s = 0.00036615946706607756	(ms)		: a0s=b0s
        zetas = 13.419130866269455		(1)
        gms = 0.14082624570054866		(1)
        smax = 5.941545585888373		(ms)
        vvh = -53.184249317587984		(mV)
        vvs = 0.7672523706054653		(mV)
        ar2=1							(1)		: 1=no inact., 0=max inact.
	ena		(mV)
	ek		(mV)
	eca		(mV)
	celsius
	v 		(mV)
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
}

ASSIGNED {
	ina 		(mA/cm2)
	ik  		(mA/cm2)
	ica 		(mA/cm2)
    ina_ina  (milliamp/cm2)  :to monitor
    ik_ik    (milliamp/cm2)  :to monitor
    ica_ica  (milliamp/cm2)  :to monitor
    itot     (milliamp/cm2)  :to monitor -- total current carried by this channel
	thegna		(mho/cm2)
	thegna_pers	(mho/cm2)
	thegna_trans	(mho/cm2)
	mpinf
	hpinf
	hptau	(ms)
	gna			(mho/cm2)
	gk			(mho/cm2)
	gca			(mho/cm2)
	erev		(mV)
	ptot
	minf
	hinf
	mtau (ms)
	htau (ms)
	sinf (ms)
	taus (ms)
}


STATE { m h s hp }

BREAKPOINT {
        SOLVE states METHOD cnexp
        : transient arm -- IDENTICAL to na12mut
        thegna_trans = gbar*m*m*m*h*s
        : persistent arm -- 0 unless pers_frac > 0
        mpinf = 1.0/(1.0 + exp((v - pers_vhalf)/(-pers_slope)))
        thegna_pers = gbar*pers_frac*mpinf*hp
        thegna = thegna_trans + thegna_pers

        ptot = pna + pk + pca
        if (ptot <= 0) { ptot = 1 }		: guard; all three g's are 0 anyway

        gna = thegna * pna/ptot
        gk  = thegna * pk /ptot
        gca = thegna * pca/ptot

        if (use_ion_erev > 0) {
            ina  = gna*(v - ena)
            ik   = gk *(v - ek)
            ica  = gca*(v - eca)
            erev = (pna*ena + pk*ek + pca*eca)/ptot
        } else {
            ina  = gna*(v - Ena)
            ik   = gk *(v - Ek)
            ica  = gca*(v - Eca)
            erev = (pna*Ena + pk*Ek + pca*Eca)/ptot
        }

        ina_ina = ina					: keeps the currentscape naming convention
        ik_ik   = ik
        ica_ica = ica
        itot    = ina + ik + ica
}

INITIAL {
	trates(v,ar2,sh)
	m=minf
	h=hinf
	s=sinf
	hp=hpinf
}


FUNCTION alpv(v) {
         alpv = 1/(1+exp((v-vvh-sh)/vvs))
}

FUNCTION alps(v) {
  alps = exp(1.e-3*zetas*(v-vhalfs-sh)*9.648e4/(8.315*(273.16+celsius)))
}

FUNCTION bets(v) {
  bets = exp(1.e-3*zetas*gms*(v-vhalfs-sh)*9.648e4/(8.315*(273.16+celsius)))
}

LOCAL mexp, hexp, sexp

DERIVATIVE states {
        trates(v,ar2,sh)
        m' = (minf-m)/mtau
        h' = (hinf-h)/htau
        s' = (sinf - s)/taus
        hp' = (hpinf - hp)/hptau
}

PROCEDURE trates(vm,a2,sh2) {
        LOCAL  a, b, c, qt
        qt=q10^((celsius-24)/10)
	a = trap0(vm,tha+sh2,Ra,qa)
	b = trap0(-vm,-tha-sh2,Rb,qa)
	mtau = 1/(a+b)/qt
        if (mtau<mmin) {
		mtau=mmin
		}
	minf = a/(a+b)

	a = trap0(vm,thi1,Rd,qd) : +sh2 raus
	b = trap0(-vm,-thi2,Rg,qg) : - sh2 raus
	htau =  1/(a+b)/qt
        if (htau<hmin) {
		htau=hmin
		}
	hinf = 1/(1+exp((vm-thinf)/qinf)): -sh2 raus
	c=alpv(vm)
        sinf = c+a2*(1-c)
        : --- persistent-arm inactivation (old Nap2 / Magistretti & Alonso) ---
        hpinf = 1.0/(1.0 + exp((vm - pers_h_vhalf)/pers_h_slope))
        hptau = (1.0/(2.88e-6*ptrap(vm + 17, 4.63)
                    + 6.94e-6*ptrap(-(vm + 64.4), 2.63)))/qt
        if (hptau < 1) { hptau = 1 }

        taus = bets(vm)/(a0s*(1+alps(vm)))
        if (taus<smax) {
		taus=smax
		}
}

FUNCTION trap0(v,th,a,q) {
	if (fabs(v-th) > 1e-6) {
	        trap0 = a * (v - th) / (1 - exp(-(v - th)/q))
	} else {
	        trap0 = a * q
 	}
}

FUNCTION ptrap(x, y) {   : vtrap from the old Napkca.mod (Nap2)
	if (fabs(x / y) < 1e-6) {
		ptrap = y * (1 - x / y / 2)
	} else {
		ptrap = x / (exp(x / y) - 1)
	}
}
