(* BNS_CTS.m 
   Wolfgang Tichy  12/2007       *)

(* compute right hand side of BSSN equations *)


(* variables *)
variables = {Psi, B[a], alphaP, Sigma, FPsi, FB[a], FalphaP ,FSigma,
              dPsi[a],   dB[a,b],   dalphaP[a],    dSigma[a],
             ddPsi[a,b],ddB[a,b,c],ddalphaP[a,b], ddSigma[a,b],
	     lPsi,lB[a],lalphaP,lSigma, FlPsi,FlB[a],FlalphaP,FlSigma,
              dlPsi[a],   dlB[a,b],   dlalphaP[a],    dlSigma[a],
             ddlPsi[a,b],ddlB[a,b,c],ddlalphaP[a,b], ddlSigma[a,b],
	     g[a,b], alpha, beta[a], K[a,b], 
             q, vRS[a], dq[a], dvRS[a,b], x, y}

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
  Cif == else,
    Cinstruction == "FirstAndSecondDerivsOf_S(box, index_lPsi, \ 
					index_dlPsi1, index_ddlPsi11);",
    Cinstruction == "FirstAndSecondDerivsOf_Sa(box, index_lB1, \
					index_dlB11, index_ddlB111);",
    Cinstruction == "FirstAndSecondDerivsOf_S(box, index_lalphaP, \
					index_dlalphaP1, index_ddlalphaP11);",
    Cinstruction == "FirstAndSecondDerivsOf_S(box, index_lSigma, \
					index_dlSigma1, index_ddlSigma11);",
  Cif == end,

  Cinstruction == "FirstDerivsOf_Sa(box, Ind(\"BNSdata_vRSx\"), \
					 Ind(\"BNSdata_vRSxx\"));",
  Cinstruction == "FirstDerivsOf_S(box,  Ind(\"BNSdata_q\"), \
			                 Ind(\"BNSdata_qx\"));",
  (* loop of all points *)
  Cinstruction == "forallpoints(box, ijk) {",

  (* Omega \times r term *)
  OmegaCrossR1 == - Omega y,
  OmegaCrossR2 == + Omega x,
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

  (* irrotational part of 3-vel in rotating frame*)
  vRI[a] == dSigma[a],

  (* vR[a] is 3-vel. in rotating frame *)
  vR[a] == vRS[a] + vRI[a],

  (* vI[a] is vel in inertial frame *)
  (* vI[a] == vR[a] + OmegaCrossR[a], *)

  (* compute square of u^0 in rotating frame *)
  oouzerosqr == alpha2 - Psi4 delta[b,c] (beta[b] + vR[b]) (beta[c] + vR[c]),
  uzerosqr == 1.0/oouzerosqr,

  (* rest mass density, pressure, and total energy density *)
  rho0 == Power[q/kappa, n],
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

  (* more terms which we need inside the stars *)
  Cif == (bi==0 || bi==3 || bi==4 || bi==5),
    dLnrho0[a] == (n/kappa) Power[q/kappa, n-1] dq[a],
    dalpha[a] == dalphaP[a]/Psi - alphaP dPsi[a]/Psi2,
    dbeta[a,b] == dB[a,b] + epsmatrix3d[b,a,3] Omega,
    dvR[a,b] == dvRS[a,b] + ddSigma[a,b],
    doouzerosqr[a] == 2 alpha dalpha[a] -
                    4 Psi3 dPsi[a] delta[b,c] *
                    (beta[b] + vR[b]) (beta[c] + vR[c]) -
                    2 Psi4 delta[b,c] *
                    (beta[b] + vR[b]) (dbeta[c,a] + dvR[c,a]),
    duzerosqr[a] == -uzerosqr^2 doouzerosqr[a],
    dLnuzerosqr[a] == duzerosqr[a]/uzerosqr,
    dLnalphaP[a] == dalphaP[a]/alphaP,
    dLnPsi[a] == dPsi[a]/Psi,
  Cif == end,

  (* decide if use non-linear ot linear equations *)
  Cif == nonlin, (* non-linear case *)

    vecLapB[a] == delta[b,c] (ddB[a,b,c] + (1/3) ddB[b,c,a]),

    (* equations for Psi, B[a], alphaP, Sigma *)
    FPsi    == delta[b,c] ddPsi[b,c] + Psi5 LBLB/(32 alpha2) +
               2Pi Psi5 rho, 
    FB[a]   == vecLapB[a] - LB[a,b] dLnalphaPsim6[b] -
               16Pi alpha Psi4 j[a],
    FalphaP == delta[b,c] ddalphaP[b,c] - alphaP (
               (7/8) Psi4 LBLB/(4 alpha2) + 2Pi Psi4 (rho+2S) ),

    Cif == (bi==0 || bi==3 || bi==4 || bi==5), (* ell. eqn. inside stars *)
      FSigma == delta[b,c] ddSigma[b,c] + 
               (vRS[a] + dSigma[a]) *
               (dLnrho0[a] + dLnuzerosqr[a]/2 + dLnalphaP[a] + 5 dLnPsi[a]), 
    Cif == else,
      FSigma  == Sigma,  (* set Sigma=0 outside stars *)
    Cif == end,

  Cif == else, (* linear case *)

    alphaP2 == alphaP*alphaP,
    alphaP3 == alphaP2*alphaP,
    Psi6    == Psi4*Psi2,
    Psi7    == Psi4*Psi3,

    gdlB == delta[a,b] dlB[a,b],
    LlB[a,b] == dlB[a,b] + dlB[b,a] -(2/3) delta[a,b] gdlB,
    LlBdo[a,b] == delta[a,c] delta[b,d] LlB[c,d], 
    LlBLlB == LlB[a,b] LlBdo[a,b],
    vecLaplB[a] == delta[b,c] (ddlB[a,b,c] + (1/3) ddlB[b,c,a]),

    (* linearized alpha == alphaP/Psi and vR[a] *)
    lvR[a] == dlSigma[a], 
    lalpha == lalphaP/Psi - alphaP lPsi/Psi2,
    (* linearized
      oouzerosqr = alpha2 - Psi4 delta[b,c] (beta[b] + vR[b]) (beta[c] + vR[c])*)
    loouzerosqr == 2alpha lalpha - 
                 4 Psi3 lPsi delta[b,c] (beta[b] + vR[b]) (beta[c] + vR[c]) -
                 Psi4 delta[b,c] 2 (lB[b] + lvR[b]) (beta[c] + vR[c]),
    luzerosqr == -uzerosqr^2 loouzerosqr,

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

    Cif == (bi==0 || bi==3 || bi==4 || bi==5), (* ell. eqn. inside stars *)
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
               (vRS[a] + dSigma[a]) * (
                lduzerosqr[a]/(2 uzerosqr) - 
                duzerosqr[a] luzerosqr /(2 uzerosqr*uzerosqr) +
                dlalphaP[a]/alphaP - dalphaP[a] lalphaP/alphaP2 +
                5 dlPsi[a]/Psi - 5 dPsi[a] lPsi/Psi2 ),
    Cif == else,
      FlSigma  == lSigma,  (* set Sigma=0 outside stars *)
    Cif == end,

  Cif == end,


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
  pr["#define Log(x)     log((double) (x)))\n"];
  pr["#define pow2(x)    ((x)*(x))\n"];
  pr["#define pow2inv(x) (1.0/((x)*(x)))\n"];
  pr["#define Cal(x,y,z) ((x)?(y):(z))\n\n"];

  pr["\n\n\n"];

  pr["void BNS_CTS(tVarList *vlFu, tVarList *vlu, \ 
		   tVarList *vlJdu, tVarList *vldu, tVarList *vlduDerivs, \
		   int nonlin)\n"];
  pr["{\n"];

  pr["double n = Getd(\"BNSdata_n\");\n"];
  pr["double kappa = Getd(\"BNSdata_kappa\");\n"];
  pr["double Omega = Getd(\"BNSdata_Omega\");\n"];
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
  prdecvarname[{vRS[a]},  "BNSdata_vRSx"];
  prdecvarname[{dq[a]},    "BNSdata_qx"];
  prdecvarname[{dvRS[a,b]},"BNSdata_vRSxx"];

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
