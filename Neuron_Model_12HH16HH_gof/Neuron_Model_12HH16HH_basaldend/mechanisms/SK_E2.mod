: SK-type calcium-activated potassium current
: Reference : Kohler et al. 1996

NEURON {
       SUFFIX SK_E2
       USEION k READ ek WRITE ik
       USEION ca READ cai
       RANGE gSK_E2bar, gSK_E2, ik, ec50, zTau
}

UNITS {
      (mV) = (millivolt)
      (mA) = (milliamp)
      (mM) = (milli/liter)
}

PARAMETER {
          v            (mV)
          gSK_E2bar = .000001 (mho/cm2)
          zTau = 1              (ms)
          : Ca sensitivity. Default 0.00043 mM = 430 nM reproduces the
          : original hard-coded value exactly, so this is INERT unless set.
          : Lowering it represents SK sitting in a VGCC nanodomain, where
          : local [Ca] exceeds the single well-mixed 0.1 um shell that
          : CaDynamics_E2 integrates. Hill coefficient 4.8 makes this a
          : switch: at the model's ~139 nM peak the channel is 0.4% active,
          : so scaling gSK_E2bar cannot make SK a player -- only ec50 can.
          ec50 = 0.00043       (mM)
          ek           (mV)
          cai          (mM)
}

ASSIGNED {
         zInf
         ik            (mA/cm2)
         gSK_E2	       (S/cm2)
}

STATE {
      z   FROM 0 TO 1
}

BREAKPOINT {
           SOLVE states METHOD cnexp
           gSK_E2  = gSK_E2bar * z
           ik   =  gSK_E2 * (v - ek)
}

DERIVATIVE states {
        rates(cai)
        z' = (zInf - z) / zTau
}

PROCEDURE rates(ca(mM)) {
          if(ca < 1e-7){
	              ca = ca + 1e-07
          }
          zInf = 1/(1 + (ec50 / ca)^4.8)
}

INITIAL {
        rates(cai)
        z = zInf
}
