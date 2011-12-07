(* BNS_set_dlnIntegEuler.m 
   Wolfgang Tichy  11/2011    *)

(* compute the updated q *)


(* variables *)
variables = {Psi, B[a], alphaP, Sigma, dSigma[a], q, wB[a], x, y,
             lnIntegEuler, dlnIntegEuler[a]}

constvariables = {OmegaCrossR[a], OmCrossR[a]}

(* compute in this order *)
tocompute = {

  Cinstruction == "FirstDerivsOf_S(box,index_BNSdata_Sigma,
                                   Ind(\"BNSdata_Sigmax\"));",
  Cif == (bi==1),
    Cinstruction == "\n
    copy_Var_at_i0_from_Box1_Box2(grid, Ind(\"BNSdata_Sigmax\"), 0,1);\n
    copy_Var_at_i0_from_Box1_Box2(grid, Ind(\"BNSdata_Sigmay\"), 0,1);\n
    copy_Var_at_i0_from_Box1_Box2(grid, Ind(\"BNSdata_Sigmaz\"), 0,1);\n
    copy_Var_at_i0_from_Box1_Box2(grid, Ind(\"BNSdata_wBx\"), 0,1);\n
    copy_Var_at_i0_from_Box1_Box2(grid, Ind(\"BNSdata_wBy\"), 0,1);\n
    copy_Var_at_i0_from_Box1_Box2(grid, Ind(\"BNSdata_wBz\"), 0,1);",
  Cif == end,
  Cif == (bi==2),
    Cinstruction == "\n
    copy_Var_at_i0_from_Box1_Box2(grid, Ind(\"BNSdata_Sigmax\"), 3,2);\n
    copy_Var_at_i0_from_Box1_Box2(grid, Ind(\"BNSdata_Sigmay\"), 3,2);\n
    copy_Var_at_i0_from_Box1_Box2(grid, Ind(\"BNSdata_Sigmaz\"), 3,2);\n
    copy_Var_at_i0_from_Box1_Box2(grid, Ind(\"BNSdata_wBx\"), 3,2);\n
    copy_Var_at_i0_from_Box1_Box2(grid, Ind(\"BNSdata_wBy\"), 3,2);\n
    copy_Var_at_i0_from_Box1_Box2(grid, Ind(\"BNSdata_wBz\"), 3,2);",
  Cif == end,

  (* loop of all points *)
  Cinstruction == "forallpoints(box, ijk) {",

  (* Omega \times r term, stays constant *)
  OmegaCrossR1 == - Omega y,
  OmegaCrossR2 == + Omega (x-xCM),
  OmegaCrossR3 == 0,

  (* Om \times r term, depends on Om,xcm args of this function *)
  OmCrossR1 == - Om y,
  OmCrossR2 == + Om (x-xcm),
  OmCrossR3 == 0,

  (* shift in rotating frame (in the inertial frame beta^i = B^i) *)
  bet[a]  == B[a] + OmCrossR[a],
  beta[a] == B[a] + OmegaCrossR[a],

  (* some abbreviations *)
  alpha  == alphaP/Psi,
  alpha2 == alpha*alpha,
  Psi2   == Psi*Psi,
  Psi4   == Psi2*Psi2,

  (**************)
  (* corotation *)
  (**************)
  Cif == ( ((bi<=1 || bi==5) && corot1) || ((bi>=2 && bi<=4) && corot2) ),

    (* vR[a] is zero for corotation *)
    (* vR[a] == 0,  ==>  Integ. Euler is h/uzero = - C  *)

    (* compute u^0 in rotating frame *)
    oouzerosqr == alpha2 - Psi4 delta[b,c] (bet[b]) (bet[c]),
    Cif == (oouzerosqr==0),
      oouzerosqr == 1,
    Cif == end,

    (* set log of integrated Euler *)
    lnIntegEuler == log[oouzerosqr],

  (****************)
  (* general case *)
  (****************)
  Cif == else,

    Psim4 == 1/Psi4,
    Psim2 == Psim4 Psi2,
    Psim6 == Psim4 Psim2,
    DSigmaUp[a] == Psim4 dSigma[a],
    w[a] == Psim6 wB[a],
    wBDown[a] == wB[a],
    wDown[a] == Psim2 wBDown[a],

    h == (n+1) q + 1,
    h2 == h h, 
    L2 == h2 + (wDown[c] + dSigma[c]) (w[c] + DSigmaUp[c]),
    uzerosqr == L2/(alpha2 h2),
    uzero == sqrt[uzerosqr],

    (* killing vec xi^i in rotating frame, xi^a = (xi^0, xi^i) = (1,0,0,0) *)
    xi[a] == 0,

    (* U and U0 from Gorgoulhon, PRD 63, 064029 *)
    (* here I compute all those using Omega, xCM *)
    U0[a] == ( beta[a] + xi[a] + w[a]/(h uzero) )/alpha,
    U0Down[a] == Psi4 U0[a],
    U[a] == DSigmaUp[a]/(alpha h uzero),
    (* Gamma factors *)
    Gamman == alpha uzero,
    Gamma0 == 1/sqrt[1 - U0[a] U0Down[a]],
    Gamma  == Gamman Gamma0 * 
              ( 1 - U0Down[a] U[a] - wDown[a] w[a]/(alpha2 h2 uzerosqr) ),

    (* set log of integrated Euler *)
    (* Om, xcm enter only the next 2 terms *)
    betxiw[a] == bet[a] + xi[a] + w[a]/(h uzero),
    betxiwDown[a] == Psi4 betxiw[a],
    lnIntegEuler == log[alpha2 - betxiw[a] betxiwDown[a] ] + 2 log[Gamma],

  Cif == end,

  Cinstruction == "} /* end of points loop */\n",

  (* compute cart derivs of lnIntegEuler *)
  (* CAUTION: this will write into the vars
     idlnIntegEuler, idlnIntegEuler+1, idlnIntegEuler+2 *)
  Cinstruction == "FirstDerivsOf_S(box, ilnIntegEuler, idlnIntegEuler);"
}


(* symmetries *)
g[a_,b_] :=  g[b,a] /; !OrderedQ[{a,b}]

(************************************************************************)
(* information for C output *)

(* information about the function 
   the C file will be in Cfunctionfile 
*)

CFunctionFile = "BNS_set_dlnIntegEuler.c"

(* the head of the function *)
BeginCFunction[] := Module[{},

  pr["#include \"sgrid.h\"\n"];
  pr["#include \"BNSdata.h\"\n"];
  pr["\n"];
  pr["#define Power(x,y) (pow((double) (x), (double) (y)))\n"];
  pr["#define Sqrt(x)    (sqrt((double) (x)))\n"];
  pr["#define Abs(x)     (fabs((double) (x)))\n"];
  pr["#define Log(x)     (log((double) (x)))\n"];
  pr["#define pow2(x)    ((x)*(x))\n"];
  pr["#define pow2inv(x) (1.0/((x)*(x)))\n"];
  pr["#define Cal(x,y,z) ((x)?(y):(z))\n\n"];

  pr["\n\n\n"];

  pr["void BNS_set_dlnIntegEuler(tGrid *grid, int ilnIntegEuler,
                                 int idlnIntegEuler, double Om, double xcm)\n"];
  pr["{\n"];

  pr["int corot1 = Getv(\"BNSdata_rotationstate1\",\"corotation\");\n"];
  pr["int corot2 = Getv(\"BNSdata_rotationstate2\",\"corotation\");\n"];
  pr["double n = Getd(\"BNSdata_n\");\n"];
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

  prdecvar[{lnIntegEuler},  "ilnIntegEuler"];
  prdecvar[{dlnIntegEuler}, "idlnIntegEuler"];

  prdecvarname[{x},      "x"];
  prdecvarname[{y},      "y"];

  (* prdecvarname[{g[a,b]},    "gxx"]; prdecvarname[{K[a,b]}, "Kxx"]; *)
  prdecvarname[{q},         "BNSdata_q"];
  prdecvarname[{Psi},       "BNSdata_Psi"];
  prdecvarname[{alphaP},    "BNSdata_alphaP"];
  prdecvarname[{B[a]},      "BNSdata_Bx"];
  prdecvarname[{wB[a]},     "BNSdata_wBx"];
  prdecvarname[{Sigma},     "BNSdata_Sigma"];
  prdecvarname[{dSigma[a]}, "BNSdata_Sigmax"];


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
