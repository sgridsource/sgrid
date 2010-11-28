(* setADMvars.m 
   Wolfgang Tichy  12/2007       *)

(* compute ADMvars from the BNSdata vars *)


(* variables *)
variables = {Psi, B[a], alphaP, Sigma, dB[a,b], dSigma[a],
	     psi, g[a,b], alpha, beta[a], K[a,b], rho, jdo[a], Sdo[a,b],
             q, wB[a], x, y}

constvariables = {OmegaCrossR[a]}

(* compute in this order *)
tocompute = {

  Cinstruction == "FirstDerivsOf_Sa(box, Ind(\"BNSdata_Bx\"), 
					 Ind(\"BNSdata_Bxx\"));",
  Cinstruction == "FirstDerivsOf_S(box, Ind(\"BNSdata_Sigma\"),
					Ind(\"BNSdata_Sigmax\"));",

  (* loop of all points *)
  Cinstruction == "forallpoints(box, ijk) {",

  (* set lapse and Psi4 *)
  alpha  == alphaP/Psi,
  alpha2 == alpha*alpha,
  Psi2   == Psi*Psi,
  Psi4   == Psi2*Psi2,

  (* Omega \times r term *)
  OmegaCrossR1 == - Omega y,
  OmegaCrossR2 == + Omega (x-xCM),
  OmegaCrossR3 == 0,

  (* shift in rotating frame *)
  beta[a] == B[a] + OmegaCrossR[a],

  (* get 1st derivs of B *)
  (* Note: if B^i = beta^i - (Omega \times r)^i 
	=> vecLapB = vecLapbeta , LB = Lbeta, 
	since the L of any Killingvec is zero *)
  gdB == delta[a,b] dB[a,b],
  LB[a,b] == dB[a,b] + dB[b,a] -(2/3) delta[a,b] gdB,
  LBdo[a,b] == delta[a,c] delta[b,d] LB[c,d], 

  (* set psi, g_ij and K_ij *)
  psi == 1, (* set ADMvars psi to one *)
  g[a,b] == Psi4 delta[a,b],
  K[a,b] == Psi4 LBdo[a,b] / (2 alpha),

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
    Psim2 == 1/Psi2,
    Psim4 == Psim2*Psim2,
    Psim6 == Psim4*Psim2,
    dSigmaUp[a] == Psim4 dSigma[a],
    w[a] == Psim6 wB[a],
    wBDown[a] == wB[a],
    wDown[a] == Psim2 wBDown[a],
    h == (n+1) q + 1,
    h2 == h*h,
    uzerosqr == (1 + (wDown[a] + dSigma[a]) (w[a] + dSigmaUp[a])/h2)/alpha2,
  Cif == end,

  (* rest mass density, pressure, and total energy density *)
  rho0 == Power[q/kappa, n],
  P    == q rho0,
  rhoE == rho0 (1 + n q),

  (* if q=0 all matter vars are zero, which can be enforced by uzerosqr=0 *)
  Cif == (q==0),
    uzerosqr == 0,
  Cif == end,

  (* set fluid vars in 3+1 *)
  rho  == alpha2 (rhoE + P) uzerosqr - P,
  jup[a] == alpha (rhoE + P) uzerosqr (vR[a]+beta[a]),

  (* set ADMvars *)
  jdo[a] == g[a,b] jup[b],
  (* Sup[a,b] == (rhoE + P) uzerosqr (vR[a]+beta[a]) (vR[b]+beta[b]) +
                  P ginv[a,b], *)
  vRplusbetado[a] == g[a,b] (vR[b]+beta[b]),
  Sdo[a,b] == (rhoE + P) uzerosqr vRplusbetado[a] vRplusbetado[b] + P g[a,b],   

  Cinstruction == "} /* end of points loop */\n"
}


(* symmetries *)
g[a_,b_]   :=  g[b,a] /; !OrderedQ[{a,b}]
K[a_,b_]   :=  K[b,a] /; !OrderedQ[{a,b}]
Sdo[a_,b_] :=  Sdo[b,a] /; !OrderedQ[{a,b}]

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

CFunctionFile = "setADMvars.c"

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

  pr["void setADMvars(tGrid *grid)\n"];
  pr["{\n"];

  pr["int corot1 = Getv(\"BNSdata_rotationstate1\",\"corotation\");\n"];
  pr["int corot2 = Getv(\"BNSdata_rotationstate2\",\"corotation\");\n"];
  pr["double n = Getd(\"BNSdata_n\");\n"];
  pr["double kappa = Getd(\"BNSdata_kappa\");\n"];
  pr["double Omega = Getd(\"BNSdata_Omega\");\n"];
  pr["double xCM = Getd(\"BNSdata_x_CM\");\n"];
  pr["\n"];

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

  prdecvarname[{Psi},     "BNSdata_Psi"];
  prdecvarname[{B[a]},    "BNSdata_Bx"];
  prdecvarname[{alphaP},  "BNSdata_alphaP"];
  prdecvarname[{Sigma},   "BNSdata_Sigma"];
  prdecvarname[{q},      "BNSdata_q"];
  prdecvarname[{wB[a]}, "BNSdata_wBx"];

  prdecvarname[{dB[a,b]}, 	"BNSdata_Bxx"];
  prdecvarname[{dSigma[a]},     "BNSdata_Sigmax"];

  prdecvarname[{x},      "x"];
  prdecvarname[{y},      "y"];

  prdecvarname[{g[a,b]}, "gxx"];
  prdecvarname[{psi}, "psi"];
  prdecvarname[{K[a,b]}, "Kxx"];
  prdecvarname[{alpha},  "alpha"];
  prdecvarname[{beta[a]},"betax"];
  prdecvarname[{rho},    "rho"];
  prdecvarname[{jdo[a]},   "jx"];
  prdecvarname[{Sdo[a,b]}, "Sxx"];

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
