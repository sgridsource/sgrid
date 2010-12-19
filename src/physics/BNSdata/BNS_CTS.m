(* BNS_CTS.m 
   Wolfgang Tichy  12/2007       *)

(* compute residuals of BNS ham, mom, alphaP and Sigma eqns *)


(* variables *)
variables = {Psi, B[a], alphaP, Sigma, FPsi, FB[a], FalphaP ,FSigma,
              dPsi[a],   dB[a,b],   dalphaP[a],    dSigma[a],
             ddPsi[a,b],ddB[a,b,c],ddalphaP[a,b], ddSigma[a,b],
	     lPsi,lB[a],lalphaP,lSigma, FlPsi,FlB[a],FlalphaP,FlSigma,
              dlPsi[a],   dlB[a,b],   dlalphaP[a],    dlSigma[a],
             ddlPsi[a,b],ddlB[a,b,c],ddlalphaP[a,b], ddlSigma[a,b],
	     g[a,b], alpha, beta[a], K[a,b], 
             q, wB[a], dq[a], dwB[a,b], VR[a], x, y, ddSigmadA2,ddlSigmadA2}

constvariables = {OmegaCrossR[a]}

(* compute in this order *)
tocompute = {

  Cif == nonlin,
    Cinstruction == "FirstAndSecondDerivsOf_S(box, index_Psi, \ 
			Ind(\"BNSdata_Psix\"), Ind(\"BNSdata_Psixx\"));",
    Cinstruction == "FirstAndSecondDerivsOf_Sa(box, index_B1, \
			Ind(\"BNSdata_Bxx\"), Ind(\"BNSdata_Bxxx\"));",
    Cinstruction == "FirstAndSecondDerivsOf_S(box, index_alphaP, \
			Ind(\"BNSdata_alphaPx\"), Ind(\"BNSdata_alphaPxx\"));",
    Cinstruction == "FirstAndSecondDerivsOf_S(box, index_Sigma, \
			Ind(\"BNSdata_Sigmax\"), Ind(\"BNSdata_Sigmaxx\"));",
    Cinstruction == "spec_Deriv2(box, 1, Sigma, ddSigmadA2);",
  Cif == else,
    Cinstruction == "FirstAndSecondDerivsOf_S(box, index_lPsi, \ 
					index_dlPsi1, index_ddlPsi11);",
    Cinstruction == "FirstAndSecondDerivsOf_Sa(box, index_lB1, \
					index_dlB11, index_ddlB111);",
    Cinstruction == "FirstAndSecondDerivsOf_S(box, index_lalphaP, \
					index_dlalphaP1, index_ddlalphaP11);",
    Cinstruction == "FirstAndSecondDerivsOf_S(box, index_lSigma, \
					index_dlSigma1, index_ddlSigma11);",
    Cinstruction == "spec_Deriv2(box, 1, lSigma, ddlSigmadA2);",
  Cif == end,

  Cinstruction == "FirstDerivsOf_Sa(box, Ind(\"BNSdata_wBx\"), \
					 Ind(\"BNSdata_wBxx\"));",
  Cinstruction == "FirstDerivsOf_S(box,  Ind(\"BNSdata_q\"), \
			                 Ind(\"BNSdata_qx\"));",
  (* loop of all points *)
  Cinstruction == "forallpoints(box, ijk) {",

  (* Omega \times r term *)
  OmegaCrossR1 == - Omega y,
  OmegaCrossR2 == + Omega (x-xCM),
  OmegaCrossR3 == 0,

  (* shift in rotating frame (in the inertial frame beta^i = B^i) *)
  beta[a] == B[a] + OmegaCrossR[a],

  (* get 1st derivs of B *)
  (* Note: if B^i = beta^i - (Omega \times r)^i 
	=> vecLapB = vecLapbeta , LB = Lbeta, 
	since the L of any Killingvec is zero *)
  gdB == delta[a,b] dB[a,b],
  LB[a,b] == dB[a,b] + dB[b,a] -(2/3) delta[a,b] gdB,
  LBdo[a,b] == delta[a,c] delta[b,d] LB[c,d], 
  LBLB == LB[a,b] LBdo[a,b],

  (* some abbreviations *)
  alpha  == alphaP/Psi,
  alpha2 == alpha*alpha,
  Psi2   == Psi*Psi,
  Psi3   == Psi*Psi2,
  Psi4   == Psi2*Psi2,
  Psi5   == Psi*Psi4,

  (* rest mass density *)
  rho0 == Power[q/kappa, n],

  (*****************************)
  (* BEGIN: corot/general case *)
  (**************)
  (* corotation *)
  (**************)
  Cif == ( ((bi<=1 || bi==5) && corot1) || ((bi>=2 && bi<=4) && corot2) ),
    (* vR[a] is zero for corotation *)
    vR[a] == 0,
    (* compute square of u^0 in rotating frame *)
    oouzerosqr == alpha2 - Psi4 delta[b,c] (beta[b] + vR[b]) (beta[c] + vR[c]),
    uzerosqr == 1.0/oouzerosqr,
  (****************)
  (* general case *)
  (****************)
  Cif == else,
    Psim1 == 1/Psi,
    Psim2 == Psim1*Psim1,
    Psim3 == Psim2*Psim1,
    Psim4 == Psim2*Psim2,
    Psim8 == Psim4*Psim4,
    Psim6 == Psi2*Psim8,
    Psim5 == Psim6*Psi,
    Psim7 == Psim8*Psi,
    Psim9 == Psim8*Psim1,
    h == (n+1) q + 1,
    h2 == h h,
    DSigmaUp[a] == Psim4 dSigma[a],
    dSigmaUp[a] == dSigma[a],
    w[a] == Psim6 wB[a],
    wBDown[a] == wB[a],
    wDown[a] == Psim2 wBDown[a],
    L2 == h2 + (wDown[c] + dSigma[c]) (w[c] + DSigmaUp[c]),
    uzerosqr == L2/(alpha2 h2),
    (* oouzerosqr == 1.0/uzerosqr, *)
    uzero == sqrt[uzerosqr],
    vR[a] == (w[a] + DSigmaUp[a])/(uzero*h) - beta[a],

    (* more terms which we need inside the stars *)
    drho0[a] == (n/kappa) Power[Abs[q/kappa], n-1] dq[a],
    dLnPsi[a] == dPsi[a]/Psi,
    dLnh[a] == (n+1) dq[a] / h,
    dLnalphaP[a] == dalphaP[a]/alphaP,
    dalpha[a] == dalphaP[a]/Psi - alphaP dPsi[a]/Psi2,
    dLnalpha[a] == dalpha[a]/alpha,
    dL2[a] == 2*(Psim8 wBDown[c] dwB[c,a] +
                 Psim6 (dwB[c,a] dSigma[c] + wB[c] ddSigma[a,c]) +
                 Psim4 dSigmaUp[c] ddSigma[a,c])  - 
              (8 Psim9 wBDown[c] wB[c] + 12 Psim7 wB[c] dSigma[c] +
               4 Psim5 dSigma[c] dSigmaUp[c]) dPsi[a] + 2 h2 dLnh[a],
    duzerosqr[a] == (dL2[a] - 2 L2 (dalpha[a]/alpha + dLnh[a]))/(alpha2 h2),
    duzero[a] == duzerosqr[a]/(2 uzero),
    dbeta[a,b] == dB[a,b] + epsmatrix3d[b,a,3] Omega,
   
    (* dLnrhozalphaPsi2oh[a] == dLnrho0[a] + dLnalphaP[a] + dLnPsi[a] - dLnh[a],
       dLnrhozalphaPsi6uz[a] == dLnrho0[a] + dLnalphaP[a] + 5 dLnPsi[a] +
                                duzero[a]/uzero, *)
    dLnalphaPsi2oh[a] == dLnalphaP[a] + dLnPsi[a] - dLnh[a],
    dLnalphaoh[a]     == dLnalphaP[a] - dLnPsi[a] - dLnh[a],
    dLnalphaPsi6uz[a] == dLnalphaP[a] + 5 dLnPsi[a] + duzero[a]/uzero,
    drho0PLUSrho0dLnalphaPsi2oh[a] == drho0[a] + rho0 dLnalphaPsi2oh[a],
    drho0PLUSrho0dLnalphaoh[a]     == drho0[a] + rho0 dLnalphaoh[a],
    drho0PLUSrho0dLnalphaPsi6uz[a] == drho0[a] + rho0 dLnalphaPsi6uz[a],

    divwB == delta [b,c] dwB[b,c],
    divbeta == delta [b,c] dbeta[b,c],
  Cif == end,
  (* END: corot/general case *)
  (***************************)

  (* V^i =: VR[a] = vR[a] *)
  VR[a] == vR[a],  (* set VR on grid equal to local vR *)

  (* rest mass density, pressure, and total energy density *)
  (* rho0 == Power[q/kappa, n], *)
  P    == q rho0,
  rhoE == rho0 (1 + n q),

  (* fluid vars in 3+1 *)
  rho  == alpha2 (rhoE + P) uzerosqr - P,
  j[a] == alpha (rhoE + P) uzerosqr (vR[a]+beta[a]),
  S    == 3P - rhoE + rho,

  (* dLnalphaPsim6[i] = \partial_i ln(alpha Psi^{-6}) 
			= \partial_i ln(alphaP Psi^{-7})
			= Psi^7 alphaP^{-1} \partial_i(alphaP Psi^{-7}) *)
  dLnalphaPsim6[a] == dalphaP[a]/alphaP - 7 dPsi[a]/Psi,

  (************************************************)
  (* decide if use non-linear or linear equations *)
  (************************************************)
  (* non-linear case: *)
  Cif == nonlin, (* non-linear case *)

    vecLapB[a] == delta[b,c] (ddB[a,b,c] + (1/3) ddB[b,c,a]),

    (* equations for Psi, B[a], alphaP, Sigma *)
    FPsi    == delta[b,c] ddPsi[b,c] + Psi5 LBLB/(32 alpha2) +
               2Pi Psi5 rho, 
    FB[a]   == vecLapB[a] - LB[a,b] dLnalphaPsim6[b] -
               16Pi alpha Psi4 j[a],
    FalphaP == delta[b,c] ddalphaP[b,c] - alphaP (
               (7/8) Psi4 LBLB/(4 alpha2) + 2Pi Psi4 (rho+2S) ),

    (*****************************)
    (* BEGIN: corot/general case *)
    (**************)
    (* corotation *)
    (**************)
    Cif == ( ((bi<=1 || bi==5) && corot1) || ((bi>=2 && bi<=4) && corot2) ),
      (* set Sigma to zero for corotation *)
      FSigma  == Sigma,
        (* Pedro's thing for inside stars (which I don't use anymore): 
           FSigma == delta[b,c] ddSigma[b,c] + 
                   (wB[a] + dSigma[a]) *
                   (dLnrho0[a] + dLnuzerosqr[a]/2 + dLnalphaP[a] + 5 dLnPsi[a]),
        *)
    (***************)
    (* genral case *)
    (***************)
    Cif == else,
      Cif == (bi==0 || bi==3 || bi==4 || bi==5), (* inside stars *)
        FSigma == rho0 delta[b,c] ddSigma[b,c] + 
                  dSigmaUp[c] drho0PLUSrho0dLnalphaPsi2oh[c] +
                  Psim2 (wB[c] drho0PLUSrho0dLnalphaoh[c] + rho0 divwB) -
                  h uzero Psi4 (rho0 divbeta +
                                beta[c] drho0PLUSrho0dLnalphaPsi6uz[c]),
      Cif == else, (* outside stars *)
        (* FSigma == delta[b,c] ddSigma[b,c], *)
        FSigma == ddSigmadA2,
      Cif == end,

    Cif == end, (* END: corot/general case *)

  (************************************************)
  (* linear case *)
  Cif == else,

    alphaP2 == alphaP*alphaP,
    alphaP3 == alphaP2*alphaP,
    Psi6    == Psi4*Psi2,
    Psi7    == Psi4*Psi3,

    gdlB == delta[a,b] dlB[a,b],
    LlB[a,b] == dlB[a,b] + dlB[b,a] -(2/3) delta[a,b] gdlB,
    LlBdo[a,b] == delta[a,c] delta[b,d] LlB[c,d], 
    LlBLlB == LlB[a,b] LlBdo[a,b],
    vecLaplB[a] == delta[b,c] (ddlB[a,b,c] + (1/3) ddlB[b,c,a]),

    (* linearized alpha == alphaP/Psi  *)
    lalpha == lalphaP/Psi - alphaP lPsi/Psi2,

    (***************************)
    (* corotation/general case *)
    (***************************)
    (**************)
    (* corotation *)
    (**************)
    Cif == ( ((bi<=1 || bi==5) && corot1) || ((bi>=2 && bi<=4) && corot2) ),
      (* Since vR[a] == 0, *)
      lvR[a] == 0, 
      (* linearized oouzerosqr, recall:
         oouzerosqr == alpha2 - Psi4 delta[b,c] (beta[b] + vR[b]) (beta[c] + vR[c])*)
      loouzerosqr == 2alpha lalpha - 
                     4 Psi3 lPsi delta[b,c] (beta[b] + vR[b]) (beta[c] + vR[c]) -
                     Psi4 delta[b,c] 2 (lB[b] + lvR[b]) (beta[c] + vR[c]),
      luzerosqr == -uzerosqr^2 loouzerosqr,
    (****************)
    (* general case *)
    (****************)
    Cif == else,
      lLnalpha == lalpha/alpha,
      (* dLnalpha[a] == dalpha[a]/alpha, *)
      dlalpha[a] == dlalphaP[a]/Psi - dalphaP[a] lPsi/Psi2 - 
                    lalphaP dPsi[a]/Psi2 - alphaP dlPsi[a]/Psi2 +
                    2 alphaP dPsi[a] lPsi/Psi3,
      ldLnalpha[a] == dlalpha[a]/alpha - lalpha dalpha[a]/alpha2,
      (* h == (n+1) q + 1, *)
      lh == 0,
      lLnh == 0,
      dlh[a] == 0,
      ldLnh[a] == 0,

      (* wB remains const under linearization *)
      lwB[a] == 0,
      dlwB[a,b] == 0,
      (* L2 == h2 + (wDown[c] + dSigma[c]) (w[c] + DSigmaUp[c]), *)
      lL2 == 2*(Psim8 wBDown[c] lwB[c] +
                   Psim6 (lwB[c] dSigma[c] + wB[c] dlSigma[c]) +
                   Psim4 dSigmaUp[c] dlSigma[c])  - 
                (8 Psim9 wBDown[c] wB[c] + 12 Psim7 wB[c] dSigma[c] +
                 4 Psim5 dSigma[c] dSigmaUp[c]) lPsi + 2 h2 lLnh,
      luzerosqr == (lL2 - 2 L2 (lalpha/alpha + lLnh))/(alpha2 h2),
      luzero == luzerosqr/(2 uzero),
      lwBDown[a] == lwB[a],
      (* dSigmaUp[a] == dSigma[a], *)
      dlSigmaUp[a] == dlSigma[a],
      Psim10 == Psim9*Psim1,
      ldL2[a] == 2*(Psim8 (lwBDown[c] dwB[c,a] + wBDown[c] dlwB[c,a]) +
                    Psim6 (dlwB[c,a] dSigma[c] + dwB[c,a] dlSigma[c] +
                           lwB[c] ddSigma[a,c] + wB[c] ddlSigma[c,a]) +
                    Psim4 (dlSigmaUp[c] ddSigma[a,c] + 
                           dSigmaUp[c] ddlSigma[a,c] )) -
                 (16 Psim9 wBDown[c] dwB[c,a] +
                  12 Psim7 (dwB[c,a] dSigma[c] + wB[c] ddSigma[a,c]) +
                  8 Psim5 dSigmaUp[c] ddSigma[a,c] ) lPsi -
                 (8 Psim9 wBDown[c] wB[c] + 12 Psim7 wB[c] dSigma[c] +
                  4 Psim5 dSigma[c] dSigmaUp[c]) dlPsi[a] +
                 (72 Psim10 wBDown[c] wB[c] + 84 Psim8 wB[c] dSigma[c] +
                  20 Psim6 dSigma[c] dSigmaUp[c]) lPsi dPsi[a] -
                 (8 Psim9 2 wBDown[c] lwB[c] +
                  12 Psim7 (lwB[c] dSigma[c] + wB[c] dlSigma[c]) +
                  8 Psim5 dlSigma[c] dSigmaUp[c]) dPsi[a] +
                 2(h2 lLnh dLnh[a] + h dlh[a]),
      lduzerosqr[a] == -2(lLnalpha + lLnh) duzerosqr[a] + 
                       (ldL2[a] - 
                        2 lL2 (dLnalpha[a] + dLnh[a]))/(alpha2 h2) - 
                       2 L2 (ldLnalpha[a] + ldLnh[a])/(alpha2 h2),
      (* vR[a] == (wB[a] + dSigma[a])/(uzero*h) - beta[a], *)
      lvR[a] == (lwB[a] + dlSigma[a])/(uzero*h) - lB[a] +
                (wB[a] + dSigma[a]) (-luzero/(uzerosqr*h) - lh/(uzero*h2)),

      divlwB == delta [b,c] dlwB[b,c],
      divlbeta == delta [b,c] dlB[b,c],
      (* ldLnalpha[a] == dlalpha[a]/alpha - lalpha dalpha[a]/alpha2, *)
      (* rho0 == Power[q/kappa, n], *)
      (* drho0[a] == (n/kappa) Power[q/kappa, n-1] dq[a], *)
      (* ldrho0[a] == (n/kappa) Power[q/kappa, n-1] dlq[a] +
                     (n(n-1)/(kappa*kappa)) Power[q/kappa, n-2] lq dq[a], *)
      lrho0 == 0,
      ldrho0[a] == 0, (* since lq == 0 *)
      ldLnalphaP[a] == dlalphaP[a]/alphaP - lalphaP dalphaP[a]/alphaP2,
      ldLnPsi[a] == dlPsi[a] Psim1 - lPsi dPsi[a] Psim2,
    Cif == end,
    (****************************)
    (* END: corot./general case *)
    (****************************)

    (* rho  == alpha2 (rhoE + P) uzerosqr - P,
       j[a] == alpha (rhoE + P) uzerosqr (vR[a]+beta[a]),
       S    == 3P - rhoE + rho, *)
    (* linearized fluid vars in 3+1 *)
    lrho  == 2 alpha lalpha (rhoE + P) uzerosqr + 
             alpha2 (rhoE + P) luzerosqr,
    lj[a] == lalpha (rhoE + P) uzerosqr (vR[a]+beta[a]) + 
             alpha (rhoE + P) luzerosqr (vR[a]+beta[a]) +
             alpha (rhoE + P) uzerosqr (lvR[a]+lB[a]),
    lS    == lrho,

    ldLnalphaPsim6[a] == dlalphaP[a]/alphaP - dalphaP[a] lalphaP/alphaP2 -
                         7 dlPsi[a]/Psi + 7 dPsi[a] lPsi/Psi2,

    (* linearized equations for Psi, B[a], alphaP, Sigma *)
    FlPsi    == delta[b,c] ddlPsi[b,c] + 7 Psi6 lPsi LBLB/(32 alphaP2) -
                (Psi7 LBLB/(16 alphaP3)) lalphaP + 
                (Psi5/(32 alpha2)) 2 LBdo[a,b] LlB[a,b] +
                2Pi 5 Psi4 lPsi rho + 2Pi Psi5 lrho, 
    FlB[a]   == vecLaplB[a] - LlB[a,b] dLnalphaPsim6[b] - 
                LB[a,b] ldLnalphaPsim6[b] - 16Pi lalphaP Psi3 j[a] -
                16Pi alphaP 3 Psi2 lPsi j[a] - 16Pi alpha Psi4 lj[a],
    FlalphaP == delta[b,c] ddlalphaP[b,c] - lalphaP (
                (-7/32) Psi6 LBLB/(alphaP2) + 2Pi Psi4 (rho+2S) ) - alphaP (
                (21/16)Psi5 lPsi LBLB/(alphaP2) +
                (7/16) (Psi4/alpha2)LBdo[a,b] LlB[a,b] + 2Pi (
                 4 Psi3 lPsi (rho+2S) + Psi4 (lrho+2lS) ) ), 

    (*****************************)
    (* BEGIN: corot/general case *)
    (**************)
    (* corotation *)
    (**************)
    Cif == ( ((bi<=1 || bi==5) && corot1) || ((bi>=2 && bi<=4) && corot2) ),
      (* set Sigma to zero for corotation *)
      FlSigma  == lSigma,
      (* Pedro's thing (which I don't use anymore): 
      /* ell. eqn. inside stars */
      Cif == (bi==0 || bi==3 || bi==4 || bi==5),
        dvR[a,b] == 0,
        dLnrho0[a] == (n/kappa) Power[q/kappa, n-1] dq[a],
        dalpha[a] == dalphaP[a]/Psi - alphaP dPsi[a]/Psi2,
        dbeta[a,b] == dB[a,b] + epsmatrix3d[b,a,3] Omega,
        doouzerosqr[a] == 2 alpha dalpha[a] -
                        4 Psi3 dPsi[a] delta[b,c] *
                        (beta[b] + vR[b]) (beta[c] + vR[c]) -
                        2 Psi4 delta[b,c] *
                        (beta[b] + vR[b]) (dbeta[c,a] + dvR[c,a]),
        duzerosqr[a] == -uzerosqr^2 doouzerosqr[a],
        dLnuzerosqr[a] == duzerosqr[a]/uzerosqr,
        dLnuzero[a] == dLnuzerosqr[a]/2,
        dLnalphaP[a] == dalphaP[a]/alphaP,
        dLnPsi[a] == dPsi[a]/Psi,

        dlalpha[a] == dlalphaP[a]/Psi - dalphaP[a] lPsi/Psi2 - 
                      lalphaP dPsi[a]/Psi2 - alphaP dlPsi[a]/Psi2 +
                      2 alphaP dPsi[a] lPsi/Psi3,
        ldoouzerosqr[a] == 2 lalpha dalpha[a] + 2 alpha dlalpha[a] -
                           (12 Psi2 lPsi dPsi[a] + 4 Psi3 dlPsi[a]) delta[b,c] *
                           (beta[b] + vR[b]) (beta[c] + vR[c]) -
                           4 Psi3 dPsi[a] delta[b,c] *
                           2 (beta[b] + vR[b]) (lB[c] + dlSigma[c]) -
                           8 Psi3 lPsi delta[b,c] *
                           (beta[b] + vR[b]) (dbeta[c,a] + dvR[c,a]) -
                           2 Psi4 delta[b,c] * (
                            (lB[b] + dlSigma[b]) (dbeta[c,a] + dvR[c,a]) + 
                            (beta[b] + vR[b]) (dlB[c,a] + ddlSigma[c,a]) ),
        lduzerosqr[a] == -2 uzerosqr luzerosqr doouzerosqr[a] - 
                          uzerosqr^2 ldoouzerosqr[a],

        FlSigma  == delta[b,c] ddlSigma[b,c] + 
                 dlSigma[a] * ( 
                  dLnrho0[a] + dLnuzerosqr[a]/2 + dLnalphaP[a] + 5 dLnPsi[a] )+
                 (wB[a] + dSigma[a]) * (
                  lduzerosqr[a]/(2 uzerosqr) - 
                  duzerosqr[a] luzerosqr /(2 uzerosqr*uzerosqr) +
                  dlalphaP[a]/alphaP - dalphaP[a] lalphaP/alphaP2 +
                  5 dlPsi[a]/Psi - 5 dPsi[a] lPsi/Psi2 ),
      Cif == else,
        FlSigma  == lSigma,
      Cif == end, 
      *)
    (****************************)
    (* genral case (not corot.) *)
    (****************************)
    Cif == else,
      Cif == (bi==0 || bi==3 || bi==4 || bi==5), (* inside stars *)
        lLnuzero == luzero/uzero,
        dLnuzero[a] == duzerosqr[a]/(2 uzerosqr),
        lduzero[a] == lduzerosqr[a]/(2 uzero) -
                      luzerosqr duzerosqr[a] / (4 uzero uzerosqr),
        ldLnuzero[a] == lduzero[a]/(uzero) - (lLnuzero) dLnuzero[a],
        lhuzeroPsi4 == lh uzero Psi4 + h luzero Psi4 + 4 h uzero Psi3 lPsi, 
        ldLnalphaPsi2oh[a] == ldLnalphaP[a] + ldLnPsi[a] - ldLnh[a],
        ldLnalphaoh[a]     == ldLnalphaP[a] - ldLnPsi[a] - ldLnh[a],
        ldLnalphaPsi6uz[a] == ldLnalphaP[a] + ldLnuzero[a] + 5 ldLnPsi[a],
        ldrho0PLUSrho0dLnalphaPsi2oh[a] == ldrho0[a] +
                                           lrho0 dLnalphaPsi2oh[a]+
                                           rho0 ldLnalphaPsi2oh[a],
        ldrho0PLUSrho0dLnalphaoh[a]     == ldrho0[a] +
                                           lrho0 dLnalphaoh[a] +
                                           rho0 ldLnalphaoh[a],
        ldrho0PLUSrho0dLnalphaPsi6uz[a] == ldrho0[a] +
                                           lrho0 dLnalphaPsi6uz[a] +
                                           rho0 ldLnalphaPsi6uz[a],
        FlSigma == rho0 delta[b,c] ddlSigma[b,c] + 
                  dlSigmaUp[c] drho0PLUSrho0dLnalphaPsi2oh[c] +
                  Psim2 (lwB[c] drho0PLUSrho0dLnalphaoh[c] + rho0 divlwB) -
                  h uzero Psi4 (rho0 divlbeta +
                                lB[c] drho0PLUSrho0dLnalphaPsi6uz[c]) +
                  lrho0 delta[b,c] ddSigma[b,c] +
                  dSigmaUp[c] ldrho0PLUSrho0dLnalphaPsi2oh[c] +
                  Psim2 (wB[c] ldrho0PLUSrho0dLnalphaoh[c] + lrho0 divwB) -
                  2 Psim3 lPsi (wB[c] drho0PLUSrho0dLnalphaoh[c] + rho0 divwB) -
                  lhuzeroPsi4 (rho0 divbeta +
                               beta[c] drho0PLUSrho0dLnalphaPsi6uz[c]) -
                  h uzero Psi4 (lrho0 divbeta + 
                                beta[c] ldrho0PLUSrho0dLnalphaPsi6uz[c]),
(*
        FSigma == rho0 delta[b,c] ddSigma[b,c] + 
                  dSigmaUp[c] drho0PLUSrho0dLnalphaPsi2oh[c] +
                  Psim2 (wB[c] drho0PLUSrho0dLnalphaoh[c] + rho0 divwB) -
                  h uzero Psi4 (rho0 divbeta +
                                beta[c] drho0PLUSrho0dLnalphaPsi6uz[c]),

FlSigma == rho0 delta[b,c] ddlSigma[b,c] + 
           dlSigmaUp[c] drho0PLUSrho0dLnalphaPsi2oh[c] -
           h luzero Psi4 (rho0 divbeta +
                          beta[c] drho0PLUSrho0dLnalphaPsi6uz[c]) -
           h uzero Psi4 beta[c] (ldLnuzero[c]),
*)
      Cif == else, (* outside stars *)
        (* FlSigma == delta[b,c] ddlSigma[b,c], *)
        FlSigma == ddlSigmadA2,
      Cif == end,

    Cif == end, (* END: corot/general case *)

  Cif == end, (* end of nonlin/linear case *)

  Cinstruction == "} /* end of points loop */\n"
}


(* symmetries *)
g[a_,b_] :=  g[b,a] /; !OrderedQ[{a,b}]
K[a_,b_] :=  K[b,a] /; !OrderedQ[{a,b}]

ginv[a_,b_] := ginv[b,a] /; !OrderedQ[{a,b}]
Kup[a_,b_]  := Kup[b,a]  /; !OrderedQ[{a,b}]

ddPsi[a_,b_]     := ddPsi[b,a]    /; !OrderedQ[{a,b}]
ddB[a_,b_,c_]    := ddB[a,c,b]    /; !OrderedQ[{b,c}]
LB[a_,b_]        := LB[b,a]       /; !OrderedQ[{a,b}]
ddalphaP[a_,b_]  := ddalphaP[b,a] /; !OrderedQ[{a,b}]
ddSigma[a_,b_]   := ddSigma[b,a]  /; !OrderedQ[{a,b}]

ddlPsi[a_,b_]     := ddlPsi[b,a]    /; !OrderedQ[{a,b}]
ddlB[a_,b_,c_]    := ddlB[a,c,b]    /; !OrderedQ[{b,c}]
LlB[a_,b_]        := LlB[b,a]       /; !OrderedQ[{a,b}]
ddlalphaP[a_,b_]  := ddlalphaP[b,a] /; !OrderedQ[{a,b}]
ddlSigma[a_,b_]   := ddlSigma[b,a]  /; !OrderedQ[{a,b}]

(************************************************************************)
(* information for C output *)

(* information about the function 
   the C file will be in Cfunctionfile 
*)

CFunctionFile = "BNS_CTS.c"

(* the head of the function *)
BeginCFunction[] := Module[{},

  pr["#include \"sgrid.h\"\n"];
  pr["#include \"BNSdata.h\"\n"];
  pr["\n"];
  pr["#define Power(x,y) (pow((double) (x), (double) (y)))\n"];
  pr["#define Abs(x)     (fabs((double) (x)))\n"];
  pr["#define Log(x)     (log((double) (x)))\n"];
  pr["#define pow2(x)    ((x)*(x))\n"];
  pr["#define pow2inv(x) (1.0/((x)*(x)))\n"];
  pr["#define Cal(x,y,z) ((x)?(y):(z))\n\n"];

  pr["\n\n\n"];

  pr["void BNS_CTS(tVarList *vlFu, tVarList *vlu, \ 
		   tVarList *vlJdu, tVarList *vldu, tVarList *vlduDerivs, \
		   int nonlin)\n"];
  pr["{\n"];

  pr["int corot1 = Getv(\"BNSdata_rotationstate1\",\"corotation\");\n"];
  pr["int corot2 = Getv(\"BNSdata_rotationstate2\",\"corotation\");\n"];
  pr["double n = Getd(\"BNSdata_n\");\n"];
  pr["double kappa = Getd(\"BNSdata_kappa\");\n"];
  pr["double Omega = Getd(\"BNSdata_Omega\");\n"];
  pr["double xCM = Getd(\"BNSdata_x_CM\");\n"];
  pr["\n"];

  pr["tGrid *grid = vlu->grid;\n"];
  pr["int bi;\n"];
  pr["\n"];

  pr["forallboxes(grid,bi)\n"];
  pr["{\n"];
  pr["tBox *box = grid->box[bi];\n"];
  pr["int ijk;\n\n"];
  pr["\n"];
  pr["\n"];
];

(* custom variable declaration
   we have to translate between the tensor names and the C variables  
   we may need special variables, for example for parameters
*)
variabledeclarations[] := Module[{},

  prdecvl[{FPsi, FB[a], FalphaP ,FSigma}, "vlFu"];
  prdecvl[{ Psi,  B[a],  alphaP,  Sigma}, "vlu"];

  prdecvl[{FlPsi, FlB[a], FlalphaP ,FlSigma}, "vlJdu"];
  prdecvl[{ lPsi, lB[a], lalphaP, lSigma}, "vldu"];
  prdecvl[{dlPsi[a],ddlPsi[a,b], dlB[a,b],ddlB[a,b,c], dlalphaP[a],ddlalphaP[a,b], dlSigma[a],ddlSigma[a,b]}, "vlduDerivs"];

  prdecvlindices[{ Psi,  B[a],  alphaP,  Sigma}, "vlu"];
  prdecvlindices[{lPsi, lB[a], lalphaP, lSigma}, "vldu"];
  prdecvlindices[{dlPsi[a],ddlPsi[a,b], dlB[a,b],ddlB[a,b,c], dlalphaP[a],ddlalphaP[a,b], dlSigma[a],ddlSigma[a,b]}, "vlduDerivs"];

  prdecvarname[{dPsi[a]},       "BNSdata_Psix"];
  prdecvarname[{ddPsi[a,b]},    "BNSdata_Psixx"];
  prdecvarname[{dB[a,b]}, 	"BNSdata_Bxx"];
  prdecvarname[{ddB[a,b,c]},    "BNSdata_Bxxx"];
  prdecvarname[{dalphaP[a]},    "BNSdata_alphaPx"];
  prdecvarname[{ddalphaP[a,b]}, "BNSdata_alphaPxx"];
  prdecvarname[{dSigma[a]},     "BNSdata_Sigmax"];
  prdecvarname[{ddSigma[a,b]},  "BNSdata_Sigmaxx"];

  prdecvarname[{x},      "x"];
  prdecvarname[{y},      "y"];

  (* prdecvarname[{g[a,b]}, "gxx"]; prdecvarname[{K[a,b]}, "Kxx"]; *)
  prdecvarname[{alpha},   "alpha"];
  prdecvarname[{beta[a]}, "betax"];
  prdecvarname[{q},       "BNSdata_q"];
  prdecvarname[{wB[a]},   "BNSdata_wBx"];
  prdecvarname[{dq[a]},   "BNSdata_qx"];
  prdecvarname[{dwB[a,b]},"BNSdata_wBxx"];
  prdecvarname[{VR[a]},   "BNSdata_VRx"];

  prdecvarname[{ddSigmadA2},   "BNSdata_SigmaXX"];
  prdecvarname[{ddlSigmadA2},  "BNSdata_lSigmaXX"];

  pr["\n"];
];    
(* auxillary variables are automatically inserted here *)

(* custom C-code which is placed directly after the variable declarations *)
InitializationCommands[] := Module[{},

  pr["/* Jetzt geht's los! */\n"];

];

(* the end or tail of the function (we need at least a }) *)
EndCFunction[] := Module[{},

  pr["} /* end of boxes */\n"];
  pr["\n\n"];
  pr["}  /* end of function */\n\n"];
];


(* to turn off optimization for development set optimizeflag = False *)
optimizeflag = True;

(* use 3d tensors the default is 3 *)
TensorEquationsDim = 3;

(************************************************************************)
(* now we are ready to go *)

<< "../../Math/MathToC/TensorEquationsToC.m"
