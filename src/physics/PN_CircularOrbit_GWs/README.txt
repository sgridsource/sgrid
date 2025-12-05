General usage comments
======================

To compile this module, add these lines to MyConfig:
----------------------------------------------------
libpaths += src/physics/ModeComputer
libpaths += src/physics/PN_CircularOrbit_GWs
libpaths += src/utility/NumericUtils

To run, start with e.g. ParFiles/q3S0.4S0.6.par
-----------------------------------------------
*see sgrid_PN_CircularOrbit_GWs.c for a brief explanation of the pars in it.
*e.g. the initial separation is not given directly, rather we specify the
 dimensionless initial orbital frequency PN_CircularOrbit_GWs_omega of the
 orbit.
*we can also give masses and spins
*PN_CircularOrbit_GWs_t1 and PN_CircularOrbit_GWs_t2 give initial and final
 time for the inspiral evolution
*if PN_CircularOrbit_GWs_t2 is chosen too large, sgrid crashes with
 "stepsize underflow in rkqs". This is nothing much to worry about! Since
 the PN formulation is breaking down when the two particles get close,
 divergencies occur, that the RK-integrator cannot handle. Simply reduce
 PN_CircularOrbit_GWs_t2 to a value before the crash and run again. Any data
 outputted before the crash should be fine.

Units
-----
*Most parfiles use units where the total mass
 M = PN_CircularOrbit_GWs_m1 + PN_CircularOrbit_GWs_m2 = 1
*If the total mass of interest differs by e.g. a factor of 3, one can simply
 scale all calculated quantities by appropriate power of this factor. E.g.
 the separation should be multiplied by 3 if M really is 3 but M=1 was used
 in the parfile.

Output
------
*In the output dir there will be the file orbit.t . Among other things it
 contains the time evolution of the separation, the spins, and the orbital
 phase Phi. If we divide Phi by 2*pi we get the number of orbits.
*We can also turn on output for GWs (see e.g. ParFiles/BNSq1S.2B2003.par):
 e.g. psi4_l2m2_s-2.t and h_l2m2_s-2.t contain the l=m=2 mode.


Specific comments about relationships to other papers:
======================================================
PN_CircularOrbit_GWs_OrbitEOMtype
---------------------------------
*the par PN_CircularOrbit_GWs_OrbitEOMtype selects the specific PN EOM we
 integrate. The different options listed in sgrid_PN_CircularOrbit_GWs.c
 correspond to EOMs from different papers. They mainly differ when the two
 particles get close, as is expected.

Comparison with arXiv:0802.1249v2
---------------------------------
My modes differ from the $h^{lm}$ in Eq. (9.3) and the 
$\hat{H}^{lm}$ in Eq. (9.4) in arXiv:0802.1249v2.

Note: Eq. (9.3) of arXiv:0802.1249v2 gives:
$$
h^{lm} = (2 \mu / D) (m \omega)^{2/3}
         \sqrt{16\pi/5)} e^{-im\psi} \hat{H}^{lm}
$$


My modes are then:
$$
h_{WT}^{lm} = h^{lm} m/[r (m \omega)^{2/3}]
            = [2 m_1 m_2 / (r D)]  \sqrt{16\pi/5)} e^{-im\psi} \hat{H}^{lm}
$$
The factor $m/[r (m \omega)^{2/3}]$ is one at leading PN order, but not in
my code.

I checked that my l=m=2 mode agrees.

BUT arXiv:0802.1249v2 has a constant l=2,m=0 mode at leading order! I don't
have that. This is a memory term that depends on the past trajectory...
